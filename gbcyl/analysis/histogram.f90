module histogram
use num_kind
use m_particledat
use xfunc_module, only: rho
implicit none
private

real(dp), save :: bin_width_
integer, save :: n_bins_
logical, save :: uniform_volume_
real(dp), allocatable, save :: boundaries_(:)
real(dp), allocatable, save :: weighted_midpoints_(:)

public :: bin_indices, bin_indices_w_init, bin_indices_uniform_width, init_histogram
public :: bin_base_area

interface bin_indices
   module procedure bin_indices_w_init, bin_indices_uniform_width
end interface

contains

subroutine init_histogram(bin_width, n_bins, uniform_volume)
  real(dp), intent(in) :: bin_width
  integer, intent(in) :: n_bins
  logical, intent(in) :: uniform_volume
  bin_width_ = bin_width
  n_bins_ = n_bins
  uniform_volume_ = uniform_volume

  if (uniform_volume) then
     allocate(boundaries_(n_bins + 1), weighted_midpoints_(n_bins))
     call boundaries_uniform_volume(bin_width_, boundaries_, &
          weighted_midpoints_)
  end if
end subroutine
  

pure function bin_base_area(i_bin)
  integer, intent(in) :: i_bin
  real(dp) :: bin_base_area
  if (uniform_volume_) then
     bin_base_area =  4 * atan(1._dp) * &
          (boundaries_(i_bin + 1)**2 - boundaries_(i_bin)**2)
  else
     bin_base_area = 4 * atan(1._dp) * &
          ((i_bin * bin_width_)**2 - ((i_bin - 1) * bin_width_)**2)
  end if
end function

!> Collects the particles having the index @p i_bin in @p bin_indices
!!
!! @param particles the array of particles to collect from
!! @param n_particles the number of particles in the array
!! @param bin_indices contains the bin indices corresponding to the
!! particles
!! @param i_bin the index of the bin to be collected
!! @param bin_particles the collected particles
!! @param n_bin_particles number of particles in the bin
!!
subroutine make_bin(particles, n_particles, bin_indices, i_bin, &
     bin_particles, n_bin_particles)
type(particledat), dimension(:), intent(in) :: particles
integer, intent(in) :: n_particles
integer, dimension(:), intent(in) :: bin_indices
integer, intent(in) :: i_bin
type(particledat), dimension(:), intent(out) :: bin_particles
integer, intent(out) :: n_bin_particles
  integer :: i
  n_bin_particles = 0
  do i = 1, n_particles
    if(bin_indices(i) == i_bin) then
      n_bin_particles = n_bin_particles + 1
      bin_particles(n_bin_particles) = particles(i)
    end if
  end do
end subroutine make_bin



subroutine bin_indices_uniform_width(particles, n_particles, xfunc, &
     bin_width, indices, offset, direction)
type(particledat), dimension(:), intent(in) :: particles
integer, intent(in) :: n_particles
real(dp) :: bin_width
real(dp), intent(in), optional :: offset
integer, intent(in), optional :: direction
interface 
  function xfunc(prtcl)
  use num_kind
  use m_particledat
  implicit none
  real(dp) :: xfunc 
  type(particledat), intent(in) :: prtcl  
  end function xfunc
end interface
integer, dimension(:), intent(out) :: indices
  integer :: i
  integer :: i_bin
  if (present(offset) .and. present(direction)) then
    do i=1, n_particles
      i_bin=int((offset + real(direction, dp) * xfunc(particles(i))) / &
           bin_width) + 1
      indices(i)=i_bin
    end do
  else 
    do i=1, n_particles
      i_bin=int(xfunc(particles(i)) / bin_width) + 1
      indices(i)=i_bin
    end do
  end if
end subroutine


!! Assumes binning_direction_ = 1!
elemental function bin_index(particle)
  type(particledat), intent(in) :: particle
  integer :: bin_index, j
  bin_index = -1
  do j = 1, size(boundaries_)
     !! Yes this is probably slow O(size(boundaries)) and e.g. a binary 
     !! search or sorting the particles could do better.
     if (boundaries_(j + 1) > sqrt(particle%x**2 + particle%y**2)) &
          then
        bin_index = j
        exit
     end if
  end do
end function


subroutine bin_indices_w_init(particles, indices)
  type(particledat), intent(in) :: particles(:)
  integer, intent(out) :: indices(size(particles))
  if (uniform_volume_) then
     if (allocated(boundaries_)) then 
        call bin_indices_uniform_volume(particles, boundaries_, indices)
     else
        stop 'histogram: bin_indices: Error: Module not initialized. Call ' // &
          'init first!'
     end if
  else
     call bin_indices(particles, size(particles), rho, bin_width_, indices) 
  end if
end subroutine


!! Assumes that boundaries are given in ascending order.
subroutine bin_indices_uniform_volume(particles, boundaries, indices)
  type(particledat), intent(in) :: particles(:)
  real(dp), intent(in) :: boundaries(:)
  integer :: indices(size(particles))
  integer :: i, j
  do i = 1, size(particles)
     do j = 1, size(boundaries)
        !! Yes this is probably slow O(size(boundaries)) and e.g. a binary 
        !! search or sorting the particles could do better.
        if (boundaries(j + 1) > sqrt(particles(i)%x**2 + particles(i)%y**2)) &
             then
           indices(i) = j
           exit
        end if
     end do
  end do
end subroutine


subroutine boundaries_uniform_volume(first_radius, boundaries, &
     weighted_midpoints)
  real(dp), intent(in) :: first_radius
  real(dp), intent(out) :: boundaries(:)
  real(dp), intent(out), optional :: weighted_midpoints(:)
  integer :: i
  i = 1
  boundaries(i) = 0._dp
  do i = 1, size(boundaries) - 1
     boundaries(i + 1) = sqrt(boundaries(i)**2 + first_radius**2)
  end do

  if (present(weighted_midpoints)) then
     call volume_weighted_midpoints(boundaries, weighted_midpoints)
  end if
end subroutine


subroutine volume_weighted_midpoints(boundaries, midpoints)
  real(dp), intent(in) :: boundaries(:)
  real(dp), intent(out) :: midpoints(:)
  integer :: i
  do i = 1, size(boundaries) - 1
     midpoints(i) = 2._dp/3._dp * (boundaries(i + 1)**3 - boundaries(i)**3) / &
          (boundaries(i + 1)**2 - boundaries(i)**2)
  end do
end subroutine

end module histogram
