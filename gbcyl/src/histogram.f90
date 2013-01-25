module histogram

contains



!! Collects the particles having the index @p i_bin in @p bin_indices 
!!
!! @p particles the array of particles to collect from
!! @p n_particles the number of particles in the array
!! @p bin_indices contains the bin indices corresponding to the particles
!! @p i_bin the index of the bin to be collected
!! @p bin_particles the collected particles
!! @p n_bin_particles number of particles in the bin
!!
subroutine make_bin(particles, n_particles, bin_indices, i_bin, bin_particles,&
n_bin_particles)
use particle
implicit none
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



subroutine bin_indices(particles, n_particles, xfunc, bin_width, indices, offset, direction)
use nrtype
use particle
implicit none
type(particledat), dimension(:), intent(in) :: particles
integer, intent(in) :: n_particles
real(dp) :: bin_width
real(dp), intent(in), optional :: offset
integer, intent(in), optional :: direction
interface 
  function xfunc(prtcl)
  use nrtype
  use particle
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
      i_bin=int((offset + real(direction, dp) * xfunc(particles(i))) / bin_width) + 1
      indices(i)=i_bin
    end do
  else 
    do i=1, n_particles
      i_bin=int(xfunc(particles(i)) / bin_width) + 1
      indices(i)=i_bin
    end do
  end if
end subroutine bin_indices

end module histogram
