!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
program analysis
  use state_reader, only: read_configuration
  use nrtype, only: dp
  use particle, only : particledat
  use orientational_ordering, only: orientation_parameter
  use psi6_module
  implicit none
  real(dp) :: radius, height
  type(particledat), dimension(:), pointer :: particles
  integer :: io_status
  integer :: n_particles
  real(dp), dimension(3) :: director
  real(dp) :: P2
  do  
    call read_configuration(particles, n_particles, radius, height, io_status)
    if (io_status < 0) then
      exit
    else
      call orientation_parameter(particles, n_particles, P2, director)
      write(*,*) psi6_bulk(particles, n_particles, 2*radius, 2*radius, height,&
        director) 
    end if
  end do
  if(associated(particles)) deallocate(particles)
end program analysis


