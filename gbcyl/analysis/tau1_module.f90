module tau1_module
use nrtype
use particle
use tau1_negative
implicit none
private

public :: tau1_routine
public :: tau1

interface tau1_routine
  subroutine max_tau1(particles, direction, tau1_max, layer_distance)
    use particle
    use nrtype
    implicit none
    type(particledat), intent(in) :: particles(:)
    real(sp), intent(in) :: direction(3)
    real(sp), intent(out) :: tau1_max
    real(sp), intent(out) :: layer_distance
  end subroutine max_tau1
end interface

interface tau1
  module procedure tau1_p
end interface

contains

real(dp) function tau1_p(particles, direction, layer_distance) result(tau1_value)
  type(particledat), intent(in) :: particles(:)
  real(sp), intent(in) :: direction(3)
  real(sp), intent(in) :: layer_distance
  real(sp) :: rs(size(particles))
  integer :: i
  do i=1, size(particles) 
    rs(i)=dot_product(real(position(particles(i)), sp), direction)
  end do
  tau1_value=tau1(rs, layer_distance) 
end function


end module
