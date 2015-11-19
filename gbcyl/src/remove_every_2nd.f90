!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
program remove_every_2nd
  use class_factory
  use num_kind, only: dp
  use particle, only : particledat
  use class_poly_box
  implicit none
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  integer, parameter :: stdin = 5
  integer, parameter :: stdout = 6
  type(factory) :: infactory
  type(factory) :: outfactory
  type(poly_box) :: simbox
  integer :: i=0
  do  
    call factory_readstate(infactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else if(mod(i, 2)==0) then
      call factory_writestate(outfactory, stdout, simbox, particles)
    end if
    i=i+1
  end do
end program


