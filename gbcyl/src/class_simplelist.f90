!> Implements a cell list for computation of short-ranged non-bonded
!! interactions between particles.
module class_simplelist
use num_kind, only: dp
use particle
use class_poly_box
implicit none
private

public :: new_simplelist
public :: simplelist_deallocate
public :: simplelist_update
public :: simplelist_nbrmask
public :: simplelist_cell_nbrmask
public :: simplelist_allocate
public :: simplelist
public :: cell_index
public :: simplelist_nbr_cells
public :: flat_index

!> Stores the cell list. 
type simplelist
   type(poly_box) :: cached_box
   real(dp) :: threshold = 0.0
   real(dp) :: min_length
   real(dp) :: min_boundary_width
   integer :: nx, ny, nz
   real(dp) :: lx, ly, lz
   logical :: is_nx_even = .false., is_ny_even = .false., is_nz_even = .false.
   integer, allocatable, dimension(:,:,:,:) :: indices
   integer, allocatable, dimension(:,:,:) :: counts
   integer, allocatable, dimension(:,:) :: coords
   real(dp), allocatable, dimension(:,:) :: xyzlist
end type simplelist

contains

!> Constructs a new cell list @p sl for @p particles in @p simbox.
!! @p min_length is the minimum cell side length. @p is_nx_even,
!! @p is_ny_even and @p is_nz_even define if an even number of cells is
!! wanted in x, y, and/or z-directions. @p threshold is the minimum
!! change in the position of a particle that can take place in either
!! the x, y, or z-direction before the cell list is updated. The update
!! threshold is computed based on @p cutoff if @p threshold is not
!! given but @p cutoff is. @p cutoff should be the cutoff radius of the
!! interactions between particles.
!!
!! Note that the resulting cell side lengths may be different in
!! different directions. 
!! 
  subroutine new_simplelist(sl, simbox, particles, min_length, &
       min_boundary_width, is_nx_even, is_ny_even, is_nz_even)
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: particles(:)
  real(dp), intent(in) :: min_length
  real(dp), intent(in) :: min_boundary_width
  logical, intent(in), optional :: is_nx_even, is_ny_even, is_nz_even
  type(simplelist), intent(out) :: sl
  sl%cached_box = simbox
  sl%min_length = min_length
  sl%min_boundary_width = min_boundary_width
  if(present(is_nx_even)) sl%is_nx_even = is_nx_even
  if(present(is_ny_even)) sl%is_ny_even = is_ny_even
  if(present(is_nz_even)) sl%is_nz_even = is_nz_even
  call calculate_dimensions(sl, simbox)
  call simplelist_allocate(sl, size(particles))
  !! :TODO: Could make the indices list size parameterizable or optimizable
  call simplelist_populate(sl, simbox, particles)
end subroutine

!> Allocates the memory for the cell list @p sl for a system containing
!! @p n particles.
subroutine simplelist_allocate(sl, n)
  type(simplelist), intent(inout) :: sl
  integer, intent(in) :: n
  allocate(sl%indices(n, 0:sl%nx - 1, 0:sl%ny - 1, 0:sl%nz - 1))
  allocate(sl%counts(0:sl%nx - 1, 0:sl%ny - 1, 0:sl%nz - 1))
  allocate(sl%coords(n, 3))
  allocate(sl%xyzlist(n, 3))
end subroutine simplelist_allocate

!> Sets the dimensions and the number of cells in @p sl.
subroutine calculate_dimensions(sl, simbox)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox

  sl%nx = max(int(getx(simbox) / sl%min_length), 1)
  sl%ny = max(int(gety(simbox) / sl%min_length), 1)
  sl%nz = max(int(getz(simbox) / sl%min_length), 1)

  if (sl%is_nx_even) sl%nx = (sl%nx / 2) * 2
  if (sl%is_ny_even) sl%ny = (sl%ny / 2) * 2
  if (sl%is_nz_even) sl%nz = (sl%nz / 2) * 2

  if (sl%nx < 3) sl%nx = 1
  if (sl%ny < 3) sl%ny = 1
  if (sl%nz < 3) sl%nz = 1

  if (sl%nx < 4 .and. isxperiodic(simbox)) sl%nx = 1
  if (sl%ny < 4 .and. isyperiodic(simbox)) sl%ny = 1
  if (sl%nz < 4 .and. iszperiodic(simbox)) sl%nz = 1

  !! Calculate cell side lengths
  sl%lx = getx(simbox) / sl%nx
  sl%ly = gety(simbox) / sl%ny
  sl%lz = getz(simbox) / sl%nz
  !! Compute threshold from cell dimensions
  sl%threshold = 0.5 * (minval([sl%lx, sl%ly, sl%lz]) - &
       sl%min_length +  sl%min_boundary_width)
end subroutine

!> Updates the cell list @p sl based on @p simbox size and the positions
!! of @p particles.
subroutine simplelist_update(sl, simbox, particles)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp) :: maxdiff
  integer :: i
  type(simplelist) :: temp
  logical :: update_needed
  maxdiff = 0.0

  !! Check if simbox should have a different decomposition compared to cached
  !! box.
  temp%min_length = sl%min_length
  call calculate_dimensions(temp, simbox)
  update_needed = temp%nx /= sl%nx .or. temp%ny /= sl%ny .or. temp%nz /= sl%nz
  
  if (.not. update_needed) then
     !! Find out the particle that has moved the most in one direction:
     !$OMP PARALLEL DO shared(sl, particles), reduction(max:maxdiff), private(i)
     do i = 1, size(particles)
        maxdiff = maxval([maxval(abs(sl%xyzlist(i, :) - &
             position(particles(i)))), maxdiff])
     end do
     !$OMP END PARALLEL DO
  end if
  
  if (update_needed .or. maxdiff > sl%threshold) then
     sl%cached_box = simbox
     call simplelist_deallocate(sl)
     call calculate_dimensions(sl, simbox)
     call simplelist_allocate(sl, size(particles))
     call simplelist_populate(sl, simbox, particles)
  end if
end subroutine

!> Frees the memory used by @p sl.
subroutine simplelist_deallocate(sl)
  type(simplelist), intent(inout) :: sl
  if (allocated(sl%indices)) deallocate(sl%indices)
  if (allocated(sl%counts)) deallocate(sl%counts)
  if (allocated(sl%coords)) deallocate(sl%coords)
  if (allocated(sl%xyzlist)) deallocate(sl%xyzlist)
end subroutine

!> Stores the cell indices of @p particles to @p sl. @p simbox
!! dimensions are used to compute the cell indices.
!!
!! @todo make this routine vectorizable.
pure subroutine simplelist_populate(sl, simbox, particles)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer :: ix, iy, iz
  integer :: iparticle
  sl%indices = 0
  sl%counts = 0
  do iparticle = 1, size(particles)
     !! Calculate cell index assuming that the coordinates are centered
     !! at the origin
     !! We also assume that the coordinates have been transformed as
     !! if (particle%x >= getx(simbox) / 2) then
     !!     particle%x = particle%x - getx(simbox)
     !! end if
     !ix = int((particles(iparticle)%x / getx(simbox) + 0.5) * sl%nx)
     !iy = int((particles(iparticle)%y / gety(simbox) + 0.5) * sl%ny)
     !iz = int((particles(iparticle)%z / getz(simbox) + 0.5) * sl%nz)
     call cell_index(sl, particles(iparticle), simbox, ix, iy, iz)
     !! add to the right position in simplelist
     sl%counts(ix, iy, iz) = sl%counts(ix, iy, iz) + 1
     sl%indices(sl%counts(ix, iy, iz), ix, iy, iz) = iparticle
     sl%coords(iparticle, :) = (/ix, iy, iz/)
     sl%xyzlist(iparticle, :) = position(particles(iparticle))
  end do
end subroutine simplelist_populate

elemental subroutine cell_index(sl, particle, simbox, ix, iy, iz)
  type(simplelist), intent(in) :: sl
  class(particledat), intent(in) :: particle
  type(poly_box), intent(in) :: simbox
  integer, intent(out) :: ix, iy, iz
  ix = int((particle%x / getx(simbox) + 0.5) * sl%nx)
  iy = int((particle%y / gety(simbox) + 0.5) * sl%ny)
  iz = int((particle%z / getz(simbox) + 0.5) * sl%nz)
end subroutine cell_index


!> Sets @p mask(j) true if @p particles(j) is a neighbour of
!! @p particles(@p i) in the @p simbox.
pure subroutine simplelist_nbrmask(sl, simbox, i, mask)
  type(simplelist), intent(in) :: sl
  type(poly_box), intent(in) :: simbox
  integer, intent(in) :: i 
  logical, intent(out) :: mask(:)
  integer :: ix,iy,iz !! cell indices of particle i
  mask=.false.
  ix = sl%coords(i,1)
  iy = sl%coords(i,2)
  iz = sl%coords(i,3)
  call simplelist_cell_nbrmask(sl, simbox, ix, iy, iz, mask)
  mask(i)=.false.
end subroutine

!> Sets @p cell_nbrmask(j) true if @p particles(j) is in cell @p ix, 
!! @p iy, @p iz or one of its nearest neighbouring cells in the cell
!! list @p sl. @p simbox is the simulation box where the particles
!! reside.
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

subroutine simplelist_nbr_cells(sl, ix, iy, iz, nbr_cells, n_nbr_cells)
  type(simplelist), intent(in) :: sl
  integer, intent(in) :: ix, iy, iz
  integer, intent(out) :: nbr_cells(3, 27)
  integer, intent(out) :: n_nbr_cells
  integer :: x, xl, y, yl, z, zl
  associate(simbox => sl%cached_box)
    n_nbr_cells = 0
    nbr_cells = -1
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
                      n_nbr_cells = n_nbr_cells + 1
                      nbr_cells(:, n_nbr_cells) = [xl, yl, zl]
                   end if
                end do
             end if
          end do
       end if
    end do
  end associate
end subroutine simplelist_nbr_cells

function flat_index(sl, ix, iy, iz) result(i)
  type(simplelist), intent(in) :: sl
  integer, intent(in) :: ix, iy, iz
  integer :: i
  ! check indices
  if (0 <= ix .and. ix < sl%nx .and. 0 <= iy .and. iy < sl%ny .and. &
       0 <= iz .and. iz < sl%nz) then
     i = ix + iy * sl%nx + iz * sl%nx * sl%ny
  else
     i = -1
  end if
end function flat_index

end module class_simplelist
