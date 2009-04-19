module tau1_module
use nrtype, only: sp, dp
use nr, only: mnbrak, brent
use orientational_ordering, only: orientation_parameter, director
use particle
implicit none



  public :: tau1_data
  public :: tau1_parameters
  public :: tau1_routine
  public :: create_parameters
  public :: calculate_and_print
  public :: tau1_f
  public :: to_string
  


PRIVATE 

  type tau1_data
    real(dp) :: value
    real(dp) :: layer_spacing
  end type tau1_data



  type tau1_parameters
    real(dp) :: min_layer_spacing
    real(dp) :: max_layer_spacing
    logical :: is_direction_set
    real(dp), dimension(3) :: direction
    real(dp) :: tolerance
  end type tau1_parameters



  interface 
    function tau1(rs, n, d)
    use nrtype, only: dp
    real(dp), dimension(:), intent(in) :: rs
    integer, intent(in) :: n
    real(dp), intent(in) :: d
    end function
  end interface


  !! This could be replaced by using optional arguments in functions. 
  interface create_parameters
    module procedure create_parameters, create_parameters_wd
  end interface create_parameters


 
contains



  function create_parameters(min_layer_spacing, max_layer_spacing, tolerance) &
  result(parameters)
  implicit none
  type(tau1_parameters) :: parameters
    real(dp), intent(in) :: min_layer_spacing
    real(dp), intent(in) :: max_layer_spacing
    real(dp), intent(in) :: tolerance 
    parameters%min_layer_spacing = min_layer_spacing
    parameters%max_layer_spacing = max_layer_spacing
    parameters%is_direction_set = .false.
    parameters%direction = 0
    parameters%tolerance = tolerance 
  end function create_parameters



  function create_parameters_wd(min_layer_spacing, max_layer_spacing, &
  direction, tolerance) &
  result(parameters)
  implicit none
  type(tau1_parameters) :: parameters
  real(dp), intent(in) :: min_layer_spacing
  real(dp), intent(in) :: max_layer_spacing
  real(dp), dimension(3) :: direction
  real(dp), intent(in) :: tolerance 
    parameters = &
      create_parameters(min_layer_spacing, max_layer_spacing, tolerance)
    parameters%direction = direction
    parameters%is_direction_set = .true.
  end function create_parameters_wd



  function to_string(a_tau1_data)
  implicit none
  character(len=78) :: to_string
  type(tau1_data), intent(in) :: a_tau1_data
     write(to_string, *) a_tau1_data
  end function to_string
  
  !! :TODO: separate string representation and printing from calculation
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
  use tau1_negative
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
    !! Here are some magic numbers
    !! :TODO: Change to adjustable
    real(dp), parameter :: interval_begin = 1.7                
    real(dp), parameter :: interval_end = interval_begin + 4.4 
    real(dp), parameter :: resolution = 1e-2*(interval_end-interval_begin) 
    real(dp) :: P2
    real(dp), dimension(:), allocatable :: rs
    real(dp), dimension(3) :: director_
    !! 1. choose direction by calculating the order parameter and the 
    !!    corresponding eigenvector
    call orientation_parameter(particles, n_particles, P2, director_)
    !! 2. find maximum value of tau1(d) 
    !! 2.1. select initial points     
    !!      (how? Must be given by user if cannot be automated robustly.)
    !! 2.2. bracket with mnbrak
    if(.not. allocated(rs)) allocate(rs(n_particles))
    do i = 1, n_particles
      position = (/particles(i)%x, particles(i)%y, particles(i)%z/)
      rs(i) = dot_product(position, director_)
    end do
    !! Calculate initial estimate for the global minimum
    tau1_max = -simple_min(interval_begin, interval_end, tau1_n, &
      resolution, layer_distance)
    ax = layer_distance - resolution
    bx = layer_distance
    cx = layer_distance + resolution
    !! Magic number.
    tol = 1e-4
    tau1_max = -brent(ax, bx, cx, tau1_n, tol, layer_distance)
    direction = director_
  end subroutine tau1_routine



  function tau1_f(particles, n_particles, parameters)
  use tau1_negative
  implicit none
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: n_particles
  type(tau1_parameters) :: parameters
  type(tau1_data) :: tau1_f
    intrinsic dot_product
    real(dp) :: ax, bx, cx
    integer :: i
    real(dp), parameter :: resolution = 1e-2
    real(dp), dimension(:), allocatable :: rs
    if(.not. parameters%is_direction_set) then
      !! Calculate director.
      parameters%direction = &
        director(orientation_parameter(particles, n_particles))
      parameters%is_direction_set = .true. 
    end if
    if(.not. allocated(rs)) allocate(rs(n_particles))
    do i = 1, n_particles
       rs(i) = dot_product(posvector(particles(i)), parameters%direction)
    end do
    !! Calculate initial estimate for the global minimum
    call init(rs, n_particles)
    tau1_f%value = -simple_min(parameters%min_layer_spacing, &
      parameters%max_layer_spacing, tau1_n, &
      resolution, tau1_f%layer_spacing)
    ax = tau1_f%layer_spacing - resolution
    bx = tau1_f%layer_spacing
    cx = tau1_f%layer_spacing + resolution
    tau1_f%value = -brent(ax, bx, cx, tau1_n, parameters%tolerance, &
      tau1_f%layer_spacing)
  end function tau1_f



  function simple_min(a, b, func, resolution, xmin)
    use tau1_negative
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
    trial_fs = (/ (func(trial_xs(i)), i = 1, n_points) /) 
    minpos = minloc(trial_fs)
    xmin = trial_xs(minpos(1))
    simple_min = minval(trial_fs)
    deallocate(trial_xs, trial_fs)
  end function simple_min



end module tau1_module
