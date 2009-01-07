module orientation_parameter_wrapper
  use particle, only: particledat
  use nrtype, only: dp
  use orientational_ordering, only: orientation_parameter
  implicit none


  contains 



  subroutine calculate_and_print(particles, n_particles, radius, height)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp), intent(in) :: radius
    real(dp), intent(in) :: height
    real(dp) :: value
    real(dp), dimension(3) :: direction
    call orientation_parameter(particles, n_particles, value, direction)
    write(*, *) value, direction
  end subroutine calculate_and_print



end module orientation_parameter_wrapper
