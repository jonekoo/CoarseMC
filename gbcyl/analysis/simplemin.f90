module simplemin
use num_kind
implicit none

contains

function simple_min(a, b, func, resolution, xmin)
  real(sp), intent(in) :: a
  real(sp), intent(in) :: b
  real(sp), intent(in) :: resolution
  real(sp), intent(out) :: xmin
  real(sp) :: simple_min
  interface
    function func(x)
      use num_kind
      implicit none
      real(sp), intent(in) :: x
      real(sp) :: func
    end function func
  end interface
  real(sp), dimension(:), allocatable :: trial_xs
  real(sp), dimension(:), allocatable :: trial_fs
  integer :: n_points 
  integer, dimension(1) :: minpos
  integer :: i
  n_points = int((b - a)/resolution) + 1
  allocate(trial_xs(n_points), trial_fs(n_points))
  trial_xs = (/ (real(i, sp), i = 0, n_points - 1) /)
  trial_xs = trial_xs*resolution + a
  trial_fs = (/ (func(trial_xs(i)), i = 1, n_points) /) 
  minpos = minloc(trial_fs)
  xmin = trial_xs(minpos(1))
  simple_min = minval(trial_fs)
  if (allocated(trial_xs)) deallocate(trial_xs)
  if (allocated(trial_fs)) deallocate(trial_fs)
end function simple_min

end module
