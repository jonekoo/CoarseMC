!! Makes a t-configuration of particles.
!! Comment: Let's not make box here because we may want different kind of
!! boxes with the t-configuration inside. 
!! 
subroutine make_tconf(particles, n_particles)
implicit none
type(particledat), dimension(:), intent(inout) :: particles
integer, intent(inout) :: n_particles
!!real(dp) :: particle_length = 4.4_dp
!!real(dp) :: particle_diameter = 1.0_dp
  !! Put one particle on the center of the box, parallel to z-axis
  particles(1) = new_particle()
  particles(2) = new_particle()
  !! Put other particle above that so that they form a t-configuration
  particles(2)%z = kappa_sigma()/2.0_dp + sigma_0()/2.0_dp
  particles(2)%ux = 1.0_dp
  particles(2)%uz = 0.0_dp
  n_particles = 2
end subroutine
