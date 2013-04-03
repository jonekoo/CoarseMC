!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
program analysis
  use class_factory
  use nrtype, only: dp
  use particle, only : particledat
  use pov
  use m_fileunit
  use class_poly_box
  implicit none
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  character(len = 7) :: ichr
  integer :: i
  type(factory) :: afactory
  type(poly_box) :: simbox
  integer, parameter :: stdin = 5
  i = 0
  do 
    call readstate(afactory, stdin, simbox, particles, io_status) 
    if (io_status < 0) then
      exit
    else
      write(ichr, '(I7.7)') i
      call povout(simbox, particles(1:size(particles)), 'povout.pov.' // ichr) 
    end if
    i = i + 1
  end do
end program analysis


