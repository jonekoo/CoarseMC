module class_simplelist
use cell
use nrtype
use particle
use class_pair_potential
use class_poly_box
use class_parameter_writer
use class_parameterizer
use all_pairs
implicit none
private

public :: new_simplelist
public :: simplelist_delete
public :: delete
!public :: scaledposition
public :: simplelist_init
public :: simplelist_writeparameters
public :: update
public :: simplelist_nbrmask
public :: simplelist_allocate
public :: simplelist
public :: simplelist_pairinteractions

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

interface simplelist_pairinteractions
  module procedure simplelist_pairinteractions, simplelist_totalpairinteractions
end interface

type simplelist
  integer :: nx, ny, nz
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

pure subroutine simplelist_populate(sl, simbox, particles)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer :: ix,iy,iz
  integer :: iparticle
  sl%indices=0
  sl%counts=0
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
end subroutine

subroutine simplelist_singleupdate(sl, simbox, particles, i)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  if (any(abs(minimage(simbox, position(particles(i))-sl%xyzlist(i,:))) > this_updatethreshold)) then
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

subroutine simplelist_totalpairinteractions(sl, simbox, particles, energy, overlap)
  type(simplelist), intent(in) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap 
  real(dp) :: singleenergy
  logical, dimension(size(particles)) :: mask
  integer :: i
  logical :: overlap_i
  if (.true.) then
    call simplelist_total_by_cell(sl, simbox, particles, energy, overlap)
  !$ else if (.true.) then !! Compiled with OpenMP
  !$ overlap = .false.
  !$OMP PARALLEL DO &
  !$OMP DEFAULT(shared) & 
  !$OMP REDUCTION(+:energy) REDUCTION(.or.:overlap) &
  !$OMP private(mask, singleenergy, overlap_i)  
  !$ do i=1, size(particles)
  !$  call simplelist_nbrmask(sl, simbox, particles, i, mask)
  !$  call maskedinteractions(mask(i:), simbox, particles(i:), 1, &
  !$       singleenergy, overlap_i)
  !$  overlap = overlap .or. overlap_i
  !$  energy = energy + singleenergy
  !$ end do
  !$OMP END PARALLEL DO
  else
    do i=1, size(particles)
      call simplelist_nbrmask(sl, simbox, particles, i, mask)
      call maskedinteractions(mask(i:), simbox, particles(i:), 1, &
      singleenergy, overlap)
      if (overlap) exit
      energy = energy + singleenergy
    end do
  end if
end subroutine


subroutine simplelist_total_by_cell(sl, simbox, particles, &
energy, overlap)
  type(simplelist), intent(in) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap 
  integer :: i, j, ix, iy, iz
  logical :: mask(size(particles))
  real(dp) :: energy_j
  integer :: helper(size(particles))
  integer, allocatable :: temp_helper(:)
  integer :: n_mask
  integer :: temp_j
  type(particledat), allocatable :: temp_particles(:)
  integer :: check(size(particles))
  helper = (/(i, i=1,size(particles))/)
  energy = 0._dp
  overlap = .false.
  !! Loop over cells. This can be thought of as looping through a 2 x 2 x 2 
  !! cube
  !! of cells.
  !$OMP PARALLEL default(shared) reduction(+:energy) reduction(.or.:overlap)&
  !$OMP& private(energy_j, i, j, mask, temp_particles, temp_j, temp_helper, n_mask)
  allocate(temp_particles(size(particles)), temp_helper(size(particles))) 
  !$OMP DO collapse(3) schedule(dynamic)
  do ix=0, sl%nx-1
  do iy=0, sl%ny-1
  do iz=0, sl%nz-1
    !!$ write(*, *) 'thread:', omp_get_thread_num(), 'calculating energy for particles in cell', ix, iy, iz
    call simplelist_nbrmask(sl, simbox, ix, iy, iz, mask)
    n_mask = count(mask)
    temp_particles(1:n_mask) = pack(particles, mask)
    temp_helper(1:n_mask) = pack(helper, mask)
    do i=1, sl%counts(ix, iy, iz)
      j = sl%indices(i, ix, iy, iz)
      ! Find position of particles(j) in temp_particles:
      do temp_j = 1, n_mask
         if(temp_helper(temp_j) == j) exit
      end do
      call pairinteractions(simbox, temp_particles(temp_j:n_mask), 1, energy_j, overlap)
      energy = energy + energy_j
      if (overlap) exit
    end do
  end do
  end do
  end do
  !$OMP END DO 
  !$OMP END PARALLEL
end subroutine

!! Mask must not include i!!
pure subroutine maskedinteractions(mask, simbox, particles, i, energy, overlap)
  logical, dimension(:), intent(in) :: mask
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  integer :: j
  real(dp) :: pairenergy
  logical :: overlap_ij
  logical :: temp_mask(size(particles))
  overlap = .false.
  energy = 0._dp
  !do j=1, size(particles)
  !  if(mask(j)) then
  !    call pairv(particles(i), particles(j), simbox, pairenergy, overlap_ij)
  !    overlap = overlap .or. overlap_ij
  !    energy = energy + pairenergy
  !  end if
  !end do
  temp_mask = .false.
  temp_mask(i) = .true.
  call pairinteractions(simbox, pack(particles, mask .or. temp_mask), i, energy, overlap)
end subroutine

!! Calculate interactions between particles and particle_i
!! 
!! @p particlei must not belong to particles!
subroutine interactions(simbox, particles, particle_i, energy, overlap)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  type(particledat), intent(in) :: particle_i
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  integer :: j
  real(dp) :: pairenergy
  logical :: overlap_ij
  real(dp) :: rij(3)
  overlap = .false.
  energy = 0._dp
  !$OMP PARALLEL
  !$OMP DO PRIVATE(overlap_ij, pairenergy) REDUCTION(+ : energy) REDUCTION(.or. : overlap)
  do j=1, size(particles)
    rij = minimage(simbox, position(particles(j)) - position(particle_i))
    call pairv(particle_i, particles(j), rij, simbox, pairenergy, overlap_ij)
    overlap = overlap .or. overlap_ij
    energy = energy + pairenergy
  end do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine




pure subroutine simplelist_pairinteractions(sl, simbox, particles, i, energy, overlap)
  type(simplelist), intent(in) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  logical :: mask(size(particles))
  call simplelist_nbrmask(sl, simbox, particles, i, mask) 
  call maskedinteractions(mask, simbox, particles, i, energy, overlap) 
end subroutine



end module
