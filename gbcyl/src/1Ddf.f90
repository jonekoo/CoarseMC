!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
!! :TODO: make maximum number of bins in histogram and the interval between
!! bins adjustable!
!!
program distribution_1D
  use nrtype, only: dp
  use particle, only : particledat
  use class_factory
  use class_poly_box
  use distribution
  use utils
  implicit none
  external g3rd
  type(particledat), dimension(:), pointer :: particles
  integer :: io_status
  integer, parameter :: maxbin = 100
  real(dp) :: delr = 0.1_dp
  real(dp), dimension(maxbin) :: histogram
  real(dp) :: bin_area
  type(factory) :: afactory
  integer, parameter :: stdin = 5
  type(poly_box) :: simbox
  real(dp), parameter :: z_axis(3)=(/0._dp, 0._dp, 1._dp/)
  do  
    call readstate(afactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else
      bin_area = volume(simbox)/getz(simbox)
      call gr1d(simbox, particles, z_axis, maxbin, bin_area, delr, histogram)  
      write(*, '(100(G16.10,1X))') histogram
    end if
  end do
  if(associated(particles)) deallocate(particles)
end program 


