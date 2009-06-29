subroutine make_sidebyside(particles, n_particles)
implicit none
type(particledat), dimension(:), intent(inout) :: particles
integer, intent(inout) :: n_particles
  particles(1) = new_particle()
  particles(2) = new_particle()
  particles(1)%x = -0.5_dp*sigma_0()
  particles(2)%x = 0.5_dp*sigma_0()
  n_particles = 2
end subroutine make_sidebyside
