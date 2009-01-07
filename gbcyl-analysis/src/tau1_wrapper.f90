module tau1_wrapper
  use particle, only: particledat
  use nrtype, only: dp
  use orientational_ordering, only: orientation_parameter
  use translational_ordering, only: tau1
  implicit none


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
    real(dp) :: o_parameter
    call orientation_parameter(particles, n_particles, o_parameter, direction)
    call tau1(particles, n_particles, direction, value, layer_distance)
    write(*, *) value, layer_distance, direction
  end subroutine calculate_and_print



end module tau1_wrapper
