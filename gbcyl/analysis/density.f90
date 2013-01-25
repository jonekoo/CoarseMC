!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
program analysis
  use class_factory
  use class_poly_box
  use nrtype, only: dp
  use particle, only : particledat
  implicit none
  type(particledat), dimension(:), pointer :: particles
  integer :: io_status
  type(poly_box) :: simbox
  integer, parameter :: stdin = 5
  type(factory) :: afactory
  do  
    call readstate(afactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else
      write(*,*) real(size(particles), dp)/volume(simbox)
    end if
  end do
  if(associated(particles)) deallocate(particles)
end program analysis


