!> Routine for computing the translational order parameter $\tau_1$ for
!! detecting a smectic(-A) phase of @p particles. @ direction sets the
!! reference direction which should be perpendicular to the layers.
!! At return, @p tau1_max and @p layer_distance contain the order
!! parameter and the average distance between the layers.
subroutine max_tau1(particles, direction, tau1_max, layer_distance)
  use particle
  use num_kind
  use simplemin
  use tau1_negative
  implicit none
  type(particledat), dimension(:), intent(in) :: particles
  real(sp), dimension(3), intent(in) :: direction
  real(sp), intent(out) :: tau1_max
  real(sp), intent(out) :: layer_distance
  integer :: n_particles
  real(sp) :: ax, cx, tol 
  integer :: i
  !! Here are some magic numbers
  !! :TODO: Change to adjustable
  real(sp), parameter :: interval_begin = 1.7_sp
  real(sp), parameter :: interval_end = interval_begin + 4.4_sp 
  real(sp), parameter :: resolution = 1e-2_sp * (interval_end-interval_begin) 
  real(sp), dimension(:), allocatable :: rs
  real(sp) :: res
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
  call init(rs)
  tau1_max = -simple_min(interval_begin, interval_end, tau1_n, resolution, &
       layer_distance)
  ax = layer_distance - resolution
  !!bx = layer_distance
  cx = layer_distance + resolution
  tol = 1e-4_sp  ! Arbitrary choice
  !!tau1_max = -brent(ax, bx, cx, tau1_n, tol, layer_distance)
  do while (.true.) 
     res = (cx - ax) / 4
     tau1_max = -simple_min(ax, cx, tau1_n, res, layer_distance)
     ax = layer_distance - res
     cx = layer_distance + res
     if (res <= tol) exit
  end do
end subroutine

