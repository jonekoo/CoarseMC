subroutine max_tau1(particles, direction, tau1_max, layer_distance)
  use particle
  use nrtype
  use nr, only: brent
  use simplemin
  use tau1_negative
  implicit none
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

