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
  use orientation_parameter_wrapper, only: calculate_and_print
  implicit none
  real(dp) :: radius, height
  type(particledat), dimension(:), pointer :: particles
  integer :: io_status
  integer :: n_particles
  do  
    call read_configuration(particles, n_particles, radius, height, io_status)
    if (io_status < 0) then
      exit
    else
      call calculate_and_print(particles, n_particles, radius, height) 
    end if
  end do
  if(associated(particles)) deallocate(particles)
end program analysis


