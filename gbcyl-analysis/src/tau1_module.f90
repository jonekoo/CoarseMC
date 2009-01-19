module tau1_module
use nrtype, only: dp
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



  subroutine tau1_routine(particles, n_particles, value, layer_distance, &
    direction)
    implicit none
    intrinsic dot_product
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp), intent(out) :: value
    real(dp), intent(out) :: layer_distance
    real(dp), dimension(3), intent(out) :: direction
    real(dp) :: ax, bx, cx, fa, fb, fc 
    integer :: i
    real(dp), dimension(3) :: position
    !! 1. choose direction by calculating the order parameter and the 
    !!    corresponding eigenvector
    call orientation_parameter(particles, n_particles, value, director_)
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
    ax = 3.9
    bx = 4.1
    call mnbrak(ax, bx, cx, fa, fb, fc, tau1_negative)
    !! (2.3. if needed, improve with brent.)
    direction = director_
    value = -fb
    layer_distance = bx
  end subroutine tau1_routine



  function tau1_negative(layer_distance) 
    implicit none
    real(dp), intent(in) :: layer_distance
    real(dp) :: tau1_negative
    tau1_negative = -tau1(rs_, n_particles_, layer_distance)   
  end function tau1_negative



end module tau1_module
