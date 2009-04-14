program povize
use io, only : povout
use state_reader
use nrtype
use particle, only : particledat
implicit none
real(dp) :: radius, height
type(particledat), dimension(:), pointer :: particles
  integer :: n_particles
  integer :: iostat
  call read_configuration(particles, n_particles, radius, height, iostat)
  call povout(particles, radius, height)
  if(associated(particles)) deallocate(particles)
end program povize
