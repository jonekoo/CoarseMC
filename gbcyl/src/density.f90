!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
program analysis
  use m_particle_factory
  use class_poly_box
  use num_kind, only: dp
  use m_particledat, only : particledat
  implicit none
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  type(poly_box) :: simbox
  integer, parameter :: stdin = 5
  type(factory) :: afactory
  do  
    call factory_readstate(afactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else
      write(*,*) real(size(particles), dp)/simbox%volume()
    end if
  end do
end program analysis


