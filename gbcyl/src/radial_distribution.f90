!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
!! :TODO: make maximum number of bins in histogram and the interval between
!! bins adjustable!
!!
program analysis
  use nrtype, only: dp
  use particle, only : particledat
  use class_factory
  use class_poly_box
  implicit none
  external g3rd
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  integer, parameter :: maxbin = 100
  real(dp) :: delr = 0.1_dp
  real(dp), dimension(maxbin) :: histogram
  type(factory) :: afactory
  integer, parameter :: stdin = 5
  type(poly_box) :: simbox
  do  
    call readstate(afactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else
      call gr3d(particles, size(particles), getx(simbox), gety(simbox), getz(simbox),&
      isxperiodic(simbox), isyperiodic(simbox), iszperiodic(simbox), maxbin, delr, histogram)  
      write(*, '(100(G16.10,1X))') histogram
    end if
  end do
end program analysis


