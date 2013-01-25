subroutine bin_indices(particles, n_particles, xfunc, bin_width, indices, offset, direction)
use nrtype
use particle
implicit none
type(particledat), dimension(:), intent(in) :: particles
integer, intent(in) :: n_particles
real(dp), intent(in) :: bin_width
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
