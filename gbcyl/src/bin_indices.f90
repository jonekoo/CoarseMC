!> Labels @p particles with bin indices.
!!
!! @param[in] particles the particles to label.
!! @param[in] xfunc the function defining the shape of the bins.
!! @param[in] bin_width the width of each bin.
!! @param[out] indices the bin index of each particle.
!! @param[in] offset the offset for binning.
!! @param[in] direction the direction of indexing.
!!
subroutine bin_indices(particles, xfunc, bin_width, indices, offset, direction)
use nrtype
use particle
implicit none
type(particledat), dimension(:), intent(in) :: particles
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
    do i = 1, size(particles)
      i_bin = int((offset + direction * xfunc(particles(i))) / bin_width) + 1
      indices(i) = i_bin
    end do
  else
    do i = 1, size(particles)
      i_bin = int(xfunc(particles(i)) / bin_width) + 1
      indices(i) = i_bin
    end do
  end if
end subroutine bin_indices
