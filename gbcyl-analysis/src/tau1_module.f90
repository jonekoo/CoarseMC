module tau1_module
use nrtype, only: sp, dp
use nr, only: mnbrak, brent
use orientational_ordering
use particle, only: particledat
implicit none



PRIVATE 

  public :: tau1_routine
  public :: calculate_and_print
  


  interface 
    function tau1(rs, n, d)
    use nrtype, only: dp
    real(dp), dimension(:), intent(in) :: rs
    integer, intent(in) :: n
    real(dp), intent(in) :: d
    end function
  end interface



  real(dp), dimension(:), allocatable :: rs_
  real(dp), dimension(3), save :: director_
  integer, save :: n_particles_
 
contains

  

  subroutine calculate_and_print(particles, n_particles, radius, height)
  implicit none
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: n_particles
  real(dp), intent(in) :: radius
  real(dp), intent(in) :: height
  real(dp) :: value
  real(dp) :: layer_distance
  real(dp), dimension(3) :: direction
    call tau1_routine(particles, n_particles, value, layer_distance, direction)
    write(*,*) value, layer_distance, direction
  end subroutine calculate_and_print



  subroutine tau1_routine(particles, n_particles, tau1_max, layer_distance, &
    direction)
    implicit none
    intrinsic dot_product
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp), intent(out) :: tau1_max
    real(dp), intent(out) :: layer_distance
    real(dp), dimension(3), intent(out) :: direction
    real(dp) :: ax, bx, cx, tol 
    integer :: i
    real(dp), dimension(3) :: position
    real(dp), parameter :: interval_begin = 1.7
    real(dp), parameter :: interval_end = interval_begin + 4.4
    real(dp), parameter :: resolution = 1e-2*(interval_end-interval_begin)
    real(dp) :: xmin
    real(dp) :: P2
    !! 1. choose direction by calculating the order parameter and the 
    !!    corresponding eigenvector
    call orientation_parameter(particles, n_particles, P2, director_)
    !! 2. find maximum value of tau1(d) 
    !! 2.1. select initial points     
    !!      (how? Must be given by user if cannot be automated robustly.)
    !! 2.2. bracket with mnbrak
    n_particles_ = n_particles
    if(.not. allocated(rs_)) allocate(rs_(n_particles_))
    do i = 1, n_particles_
      position = (/particles(i)%x, particles(i)%y, particles(i)%z/)
      rs_(i) = dot_product(position, director_)
    end do
    !! Calculate initial estimate for the global minimum
    tau1_max = -simple_min(interval_begin, interval_end, tau1_negative, &
      resolution, layer_distance)
    bx = layer_distance
    ax = layer_distance - resolution
    tol = 1e-4
    tau1_max = -brent(ax, bx, cx, tau1_negative, tol, layer_distance)
    direction = director_
  end subroutine tau1_routine



  function tau1_negative(layer_distance) 
    implicit none
    real(dp), intent(in) :: layer_distance
    real(dp) :: tau1_negative
    tau1_negative = -tau1(rs_, n_particles_, layer_distance)   
  end function tau1_negative



  function simple_min(a, b, func, resolution, xmin)
    use nrtype, only : sp
    implicit none
    real(sp), intent(in) :: a
    real(sp), intent(in) :: b
    real(sp), intent(in) :: resolution
    real(sp), intent(out) :: xmin
    real(sp) :: simple_min
    interface
      function func(x)
        use nrtype
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
    trial_xs = (/ (i, i = 0, n_points - 1) /)
    trial_xs = trial_xs*resolution + a
    trial_fs = (/ (tau1_negative(trial_xs(i)), i = 1, n_points) /) 
    minpos = minloc(trial_fs)
    xmin = trial_xs(minpos(1))
    simple_min = minval(trial_fs)
    deallocate(trial_xs, trial_fs)
  end function simple_min



end module tau1_module
