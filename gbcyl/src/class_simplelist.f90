module class_simplelist
use cell
use nrtype
use particle
use class_poly_box
use class_parameter_writer
use class_parameterizer
implicit none
private

public :: new_simplelist
public :: simplelist_delete
public :: delete
!public :: scaledposition
public :: simplelist_init
public :: simplelist_writeparameters
public :: update
public :: nbrmask
public :: simplelist
public :: simplelist_allocate

real(dp), save :: this_minlength = 7.5_dp
logical, save :: this_iseven = .false.
real(dp), save :: this_updatethreshold = 0.5_dp

interface simplelist_init
  module procedure initwtparameters, initwtreader 
end interface

interface new_celllist
  module procedure new_simplelist
end interface

interface update
  module procedure simplelist_update, simplelist_singleupdate !, simplelist_updatei
end interface

type simplelist
  integer :: nx, ny, nz
  integer, pointer, dimension(:,:,:,:) :: indices
  integer, pointer, dimension(:,:,:) :: counts
  integer, pointer, dimension(:,:) :: coords
  real(dp), pointer, dimension(:,:) :: xyzlist
end type

interface nbrmask
  module procedure simplelist_nbrmask
end interface

interface delete
  module procedure simplelist_delete
end interface

contains

subroutine initwtreader(parameterreader)
  type(parameterizer), intent(in) :: parameterreader
  call getparameter(parameterreader, 'cellminlength', this_minlength)
  call getparameter(parameterreader, 'isdivisioneven', this_iseven)
  call getparameter(parameterreader, 'updatethreshold', this_updatethreshold)
end subroutine

subroutine initwtparameters(minlength, iseven, updatethreshold)
  real(dp), intent(in) :: minlength
  logical, intent(in) :: iseven
  real(dp), intent(in) :: updatethreshold
  this_minlength = minlength
  this_iseven = iseven
  this_updatethreshold = updatethreshold
end subroutine

subroutine simplelist_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call writecomment(writer, 'simplelist parameters')
  call writeparameter(writer, 'cellminlength', this_minlength)
  call writeparameter(writer, 'isdivisioneven', this_iseven) 
  call writeparameter(writer, 'updatethreshold', this_updatethreshold)
end subroutine

!! Constructs a new cell list of @p particles. Given simulation box 
!! dimensions in @p simbox and the minimum cell side length in @p minlength
!! the new list is returned in the variable cl. Note that cell side lengths 
!! may be different in different directions. 
!! 
!! @p cl is the cell list. 
!! @p simbox is the simulation box. 
!! @p particles are the particles to be distributed to the cells.
!! @p nparticles is the number of particles.
!! @p minlength is the minimum length of a cell side. 
!!
function new_simplelist(simbox, particles) result(sl)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  type(simplelist) :: sl
  call calculate_dimensions(simbox, sl%nx, sl%ny, sl%nz)
  call simplelist_allocate(sl, size(particles))
  !! :TODO: Make the indices list size parameterizable or optimizable
  call simplelist_populate(sl, simbox, particles)
end function

subroutine simplelist_allocate(sl, n)
  type(simplelist), intent(inout) :: sl
  integer, intent(in) :: n
  allocate(sl%indices(n, 0:sl%nx-1, 0:sl%ny-1, 0:sl%nz-1))
  allocate(sl%counts(0:sl%nx-1,0:sl%ny-1,0:sl%nz-1))
  allocate(sl%coords(n,3))
  allocate(sl%xyzlist(n,3))
end subroutine

subroutine calculate_dimensions(simbox, nx, ny, nz)
  type(poly_box), intent(in) :: simbox
  integer, intent(out) :: nx, ny, nz
  nx = max(ncells(getx(simbox), this_minlength), 1)
  ny = max(ncells(gety(simbox), this_minlength), 1)
  nz = max(ncells(getz(simbox), this_minlength), 1)
  if (this_iseven) then
    if (nx < 2 .or. ny < 2 .or. nz < 2) then
      stop 'Could not create a cell list with even number of cells in'//& 
      &' all directions'
    end if
    nx = (nx / 2) * 2
    ny = (ny / 2) * 2
    nz = (nz / 2) * 2
  end if 
  if (nx < 3) nx=1
  if (ny < 3) ny=1
  if (nz < 3) nz=1
  if (nx < 4 .and. isxperiodic(simbox)) nx=1
  if (ny < 4 .and. isyperiodic(simbox)) ny=1
  if (nz < 4 .and. iszperiodic(simbox)) nz=1
end subroutine

subroutine simplelist_update(sl, simbox, particles)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer :: nx, ny, nz
  call calculate_dimensions(simbox, nx, ny, nz)
  if (nx/=sl%nx .or. ny/=sl%ny .or. nz/=sl%nz) then
    call delete(sl)
    sl%nx=nx
    sl%ny=ny
    sl%nz=nz
    call simplelist_allocate(sl, size(particles))
  end if
  call simplelist_populate(sl, simbox, particles)
end subroutine

subroutine simplelist_delete(sl)
  type(simplelist), intent(inout) :: sl
  deallocate(sl%indices)
  deallocate(sl%counts)
  deallocate(sl%coords)
  deallocate(sl%xyzlist)
end subroutine

subroutine simplelist_populate(sl, simbox, particles)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer :: ix,iy,iz
  integer :: iparticle
  sl%indices=0
  sl%counts=0
  !call delete(sl) 
  !sl=new_simplelist(simbox, particles)
  do iparticle = 1, size(particles)
    !! calculate cell index assuming centered coordinates
    ix=int(((particles(iparticle)%x/getx(simbox)+0.5_dp)*sl%nx))
    iy=int(((particles(iparticle)%y/gety(simbox)+0.5_dp)*sl%ny))
    iz=int(((particles(iparticle)%z/getz(simbox)+0.5_dp)*sl%nz))
    !! add to the right position in simplelist
    sl%counts(ix,iy,iz)=sl%counts(ix,iy,iz)+1
    sl%indices(sl%counts(ix,iy,iz), ix, iy, iz)=iparticle
    sl%coords(iparticle,:)=(/ix, iy, iz/)
    sl%xyzlist(iparticle,:)=position(particles(iparticle))
  end do

  !! Let's check that all the particles are in the list, although it is quite 
  !! unlikely to happen when the list is used correctly.
  if(size(particles)/=sum(sl%counts)) stop 'Some particles are not in the cell list' 
end subroutine

subroutine simplelist_singleupdate(sl, simbox, particles, i)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  if (any(abs(minimage(simbox, position(particles(i))-sl%xyzlist(i,:))) > this_updatethreshold)) then
    !write(*, *) 'called update from singleupdate'
    !write(*, *) abs(minimage(simbox, position(particles(i))-sl%xyzlist(i,:))) > this_updatethreshold
    call update(sl, simbox, particles)
  end if
end subroutine

function maskarray(n, periodic, i)
  integer, intent(in) :: n
  logical, intent(in) :: periodic 
  integer, intent(in) :: i
  logical, dimension(0:n-1) :: maskarray
  integer :: j
  if (periodic) then
    do j=i+2, i+n-2
       maskarray(mod(i,n))=.false.
    end do
  else
    maskarray(:i-2)=.false.
    maskarray(i+2:)=.false.
  end if
end function

function maskcells(nx, ny, nz, xperiodic, yperiodic, zperiodic, ix, iy, iz)
  integer, intent(in) :: nx, ny, nz
  logical, intent(in) :: xperiodic, yperiodic, zperiodic
  integer, intent(in) :: ix, iy, iz
  logical, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: maskcells
  integer :: x, y, z
  maskcells=.true.
  
  if (xperiodic) then
    do x=ix+2, ix+nx-2
       maskcells(mod(ix,nx), :, :)=.false.
    end do
  else
    maskcells(:ix-2,:,:)=.false.
    maskcells(ix+2:,:,:)=.false.
  end if

  if (yperiodic) then
    do y=iy+2, iy+ny-2
       maskcells(:, mod(iy,ny), :)=.false.
    end do
  else
    maskcells(:,:iy-2,:)=.false.
    maskcells(:,iy+2:,:)=.false.
  end if

  if (zperiodic) then
    do z=iz+2, iz+nz-2
       maskcells(:, :, mod(iz,nz))=.false.
    end do
  else
    maskcells(:,:,:iz-2)=.false.
    maskcells(:,:,iz+2:)=.false.
  end if

end function

function simplelist_nbrmask(sl, simbox, particles, i)
  type(simplelist), intent(in) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i 
  integer :: ix,iy,iz !! cell indices of particle i
  integer :: x,y,z    !! iteration indices through neighbouring cell
  integer :: xl,yl,zl
  logical, dimension(size(particles)) :: simplelist_nbrmask
  simplelist_nbrmask=.false.
  ix=sl%coords(i,1)
  iy=sl%coords(i,2)
  iz=sl%coords(i,3)
  do x=ix-1,ix+1
    xl=x
    if (isxperiodic(simbox)) xl=mod(sl%nx+xl,sl%nx)
    if (xl>=0 .and. xl<=sl%nx-1) then
      do y=iy-1,iy+1
        yl=y
        if (isyperiodic(simbox)) yl=mod(sl%ny+yl,sl%ny)
        if (yl>=0 .and. yl<=sl%ny-1) then
          do z=iz-1,iz+1
            zl=z
            if (iszperiodic(simbox)) zl=mod(sl%nz+zl,sl%nz)
            if (zl>=0 .and. zl<=sl%nz-1) then
              simplelist_nbrmask(sl%indices(1:sl%counts(xl,yl,zl),xl,yl,zl))=.true.
            end if
          end do
        end if
      end do
    end if
  end do
  simplelist_nbrmask(i)=.false.
end function 

!pure function scaledposition(simbox, particle) result(s)
!  type(poly_box), intent(in) :: simbox
!  type(particledat), intent(in) :: particle
!  real(dp), dimension(3) :: s
!  real(dp), dimension(3) :: r
!  r = position(particle)
!  s = (/r(1) / getx(simbox), r(2) / gety(simbox), r(3) / getz(simbox)/)
!end function

end module
