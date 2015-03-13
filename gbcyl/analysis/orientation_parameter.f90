!> Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
program analysis
  use class_factory
  use nrtype, only: dp
  use particle, only : particledat
  use class_poly_box
  use utils
  use orientational_ordering
  implicit none
  !real(dp) :: radius, height
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  integer, parameter :: stdin = 5
  type(factory) :: afactory
  real(dp) :: p2
  real(dp), dimension(3) :: director
  type(poly_box) :: simbox
  do  
    call factory_readstate(afactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else
      call orientation_parameter(pack(particles, particles%rod), count(particles%rod), p2, director)
      write(*, '( 4('// fmt_char_dp() //',1X))') p2, director
    end if
  end do
end program analysis


