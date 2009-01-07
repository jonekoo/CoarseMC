module cylindrical_orientation_wrapper
  use particle, only: particledat
  use nrtype, only: dp
  use cylindrical_orientation, only: eigens
  implicit none


  contains 


  !! Calculates and prints the eigenvalues and -vectors of the cylindrical
  !! orientation tensor of @p particles.
  !!
  !! @p particles the particles over which the cylindrical orientation tensor 
  !! is evalued.
  !! @p n_particles the number of particles
  !! @p radius for interface conformance only. Not used.
  !! @p height for interface conformance only. Not used.
  !!
  subroutine calculate_and_print(particles, n_particles, radius, height)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp), intent(in) :: radius
    real(dp), intent(in) :: height
    real(dp), dimension(3) :: values
    real(dp), dimension(3, 3) :: vectors
    call eigens(particles, n_particles, values, vectors)
    write(*, *) values(1:3), vectors(1:3, 1), vectors(1:3, 2), vectors(1:3, 3) 
  end subroutine calculate_and_print



end module cylindrical_orientation_wrapper
