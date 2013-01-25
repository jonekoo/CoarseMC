module tau1_module
use nrtype, only: sp
use nr, only: brent
use particle
use tau1_negative
use simplemin
implicit none
private

public :: tau1_routine
public :: tau1

!!interface 
!!  function tau1(rs, n, d)
!!  use nrtype, only: sp
!!  real(sp), dimension(:), intent(in) :: rs
!!  integer, intent(in) :: n
!!  real(sp), intent(in) :: d
!!  end function
!!end interface

interface tau1_routine
  module procedure max_tau1
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
  tau1_value=tau1(rs, size(particles), layer_distance) 
end function

subroutine max_tau1(particles, direction, tau1_max, layer_distance)
  type(particledat), dimension(:), intent(in) :: particles
  real(sp), dimension(3), intent(in) :: direction
  real(sp), intent(out) :: tau1_max
  real(sp), intent(out) :: layer_distance
  integer :: n_particles
  real(sp) :: ax, bx, cx, tol 
  integer :: i
  !! Here are some magic numbers
  !! :TODO: Change to adjustable
  real(sp), parameter :: interval_begin = 1.7_sp
  real(sp), parameter :: interval_end = interval_begin + 4.4_sp 
  real(sp), parameter :: resolution = 1e-2_sp * (interval_end-interval_begin) 
  real(sp), dimension(:), allocatable :: rs
  n_particles = size(particles)
  !! 2. find maximum value of tau1(d) 
  !! 2.1. select initial points     
  !!      (how? Must be given by user if cannot be automated robustly.)
  !! 2.2. bracket with mnbrak
  if(.not. allocated(rs)) allocate(rs(n_particles))
  do i = 1, n_particles
    rs(i) = dot_product(real(position(particles(i)),sp), direction)
  end do
  !! Calculate initial estimate for the global minimum
  call init(rs, n_particles)
  tau1_max = -simple_min(interval_begin, interval_end, tau1_n, resolution, &
  layer_distance)
  ax = layer_distance - resolution
  bx = layer_distance
  cx = layer_distance + resolution
  !! Magic number.
  tol = 1e-4_sp
  tau1_max = -brent(ax, bx, cx, tau1_n, tol, layer_distance)
end subroutine

end module
