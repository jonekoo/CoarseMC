module tau1_module
use nrtype, only: dp
use nr, only: mnbrak, brent
use orientational_ordering
use particle, only: particledat
implicit none

PRIVATE 

  public :: tau1

  real(dp), dimension(3), save :: director_
  type(particledat), dimension(:), save :: particles_
  integer, save :: n_particles_

contains



  subroutine tau1(particles, n_particles, value, layer_distance, direction)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp), intent(out) :: value
    real(dp), intent(out) :: layer_distance
    real(dp), intent(out) :: direction
    real(dp) :: ax, bx, cx, fa, fb, fc 
    particles_ = particles
    n_particles_ = n_particles
    !! 1. choose direction by calculating the order parameter and the 
    !!    corresponding eigenvector
    call orientation_parameter(particles, n_particles, value, director_)
    !! 2. find maximum value of tau1(d) 
    !! 2.1. select initial points     
    !!      (how? Must be given by user if cannot be automated robustly.)
    !! 2.2. bracket with mnbrak
    ax = 1.0
    bx = 2.0
    call mnbrak(ax, bx, cx, fa, fb, fc, -tau1)
    !! (2.3. if needed, improve with brent.)
    direction = director_
    value = fb
    layer_distance = bx
  end subroutine tau1



  function tau1(d) result(value)
    implicit none
    
  end function

end module tau1_module
