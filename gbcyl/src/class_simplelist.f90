module class_simplelist
use nrtype
use particle
use class_poly_box
implicit none
private

public :: new_simplelist
public :: simplelist_delete
public :: delete
public :: update
public :: simplelist_nbrmask
public :: simplelist_allocate
public :: simplelist

interface new_celllist
  module procedure new_simplelist
end interface

interface update
  module procedure simplelist_update
end interface

type simplelist
  real(dp) :: min_length
  integer :: nx, ny, nz
  logical :: is_x_even = .false., is_y_even = .false., is_z_even = .false.
  integer, allocatable, dimension(:,:,:,:) :: indices
  integer, allocatable, dimension(:,:,:) :: counts
  integer, allocatable, dimension(:,:) :: coords
  real(dp), allocatable, dimension(:,:) :: xyzlist
end type

interface simplelist_nbrmask
  module procedure simplelist_nbrmask, simplelist_cell_nbrmask
end interface

interface delete
  module procedure simplelist_delete
end interface

contains

!! Constructs a new cell list of @p particles. Given simulation box 
!! dimensions in @p simbox and the minimum cell side length in @p minlength
!! the new list is returned in the variable cl. Note that cell side lengths 
!! may be different in different directions. 
!! 
!! @p cl is the cell list. 
!! @p simbox is the simulation box. 
!! @p particles are the particles to be distributed to the cells.
!! @p nparticles is the number of particles.
!! @p min_length is the minimum length of a cell side. 
!!
function new_simplelist(simbox, particles, min_length, is_x_even, is_y_even, is_z_even) result(sl)
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: particles(:)
  real(dp), intent(in) :: min_length
  logical, intent(in), optional :: is_x_even, is_y_even, is_z_even
  type(simplelist) :: sl
  sl%min_length = min_length
  if(present(is_x_even)) sl%is_x_even = is_x_even
  if(present(is_y_even)) sl%is_y_even = is_y_even
  if(present(is_z_even)) sl%is_z_even = is_z_even
  call calculate_dimensions(sl, simbox)
  call simplelist_allocate(sl, size(particles))
  !! :TODO: Could make the indices list size parameterizable or optimizable
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

subroutine calculate_dimensions(sl, simbox)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox

  sl%nx = max(int(getx(simbox)/sl%min_length), 1)
  sl%ny = max(int(gety(simbox)/sl%min_length), 1)
  sl%nz = max(int(getz(simbox)/sl%min_length), 1)

  if (sl%is_x_even) sl%nx = (sl%nx / 2) * 2
  if (sl%is_y_even) sl%ny = (sl%ny / 2) * 2
  if (sl%is_z_even) sl%nz = (sl%nz / 2) * 2

  if (sl%nx < 3) sl%nx = 1
  if (sl%ny < 3) sl%ny = 1
  if (sl%nz < 3) sl%nz = 1

  if (sl%nx < 4 .and. isxperiodic(simbox)) sl%nx = 1
  if (sl%ny < 4 .and. isyperiodic(simbox)) sl%ny = 1
  if (sl%nz < 4 .and. iszperiodic(simbox)) sl%nz = 1
end subroutine

subroutine simplelist_update(sl, simbox, particles)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  call calculate_dimensions(sl, simbox)
  call simplelist_delete(sl)
  call simplelist_allocate(sl, size(particles))
  call simplelist_populate(sl, simbox, particles)
end subroutine

subroutine simplelist_delete(sl)
  type(simplelist), intent(inout) :: sl
  deallocate(sl%indices)
  deallocate(sl%counts)
  deallocate(sl%coords)
  deallocate(sl%xyzlist)
end subroutine

pure subroutine simplelist_populate(sl, simbox, particles)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer :: ix,iy,iz
  integer :: iparticle
  sl%indices = 0
  sl%counts = 0
  do iparticle = 1, size(particles)
    !! calculate cell index assuming centered coordinates
    ix = int(((particles(iparticle)%x / getx(simbox) + 0.5_dp) * &
      real(sl%nx, dp)))
    iy = int(((particles(iparticle)%y / gety(simbox) + 0.5_dp) * &
      real(sl%ny, dp)))
    iz = int(((particles(iparticle)%z / getz(simbox) + 0.5_dp) * & 
      real(sl%nz, dp)))
    !! add to the right position in simplelist
    sl%counts(ix,iy,iz)=sl%counts(ix,iy,iz)+1
    sl%indices(sl%counts(ix,iy,iz), ix, iy, iz)=iparticle
    sl%coords(iparticle,:)=(/ix, iy, iz/)
    sl%xyzlist(iparticle,:)=position(particles(iparticle))
  end do
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

pure subroutine simplelist_nbrmask(sl, simbox, particles, i, mask)
  type(simplelist), intent(in) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i 
  logical, intent(out) :: mask(:)
  integer :: ix,iy,iz !! cell indices of particle i
  mask=.false.
  ix=sl%coords(i,1)
  iy=sl%coords(i,2)
  iz=sl%coords(i,3)
  call simplelist_cell_nbrmask(sl, simbox, ix, iy, iz, mask)
  mask(i)=.false.
end subroutine

pure subroutine simplelist_cell_nbrmask(sl, simbox, ix, iy, iz, cell_nbrmask)
  implicit none
  type(simplelist), intent(in) :: sl
  type(poly_box), intent(in) :: simbox
  integer, intent(in) :: ix,iy,iz !! cell indices
  logical, intent(out) :: cell_nbrmask(:)
  integer :: x, y, z    !! iteration indices through neighbouring cell
  integer :: xl, yl, zl
  cell_nbrmask=.false.
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
              cell_nbrmask(sl%indices(1:sl%counts(xl,yl,zl),xl,yl,zl))=.true.
            end if
          end do
        end if
      end do
    end if
  end do
end subroutine

end module
