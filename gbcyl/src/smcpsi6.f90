!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Calculates and prints the bond-orientational ordering parameter psi6 for a
!! Sm-C or Sm-A phase. Other information printed is the layer distance, the 
!! orientational ordering parameter of the local layer normalss and the 
!! reference direction, the global layer normal used in the calculation of 
!! psi6.
!!
!! Usage: smctau1 < configurationfile
!!
program smcpsi6
use particle
use class_poly_box
use utils
use nrtype
!use class_cylformatter
use class_factory
use layernormal
use class_poly_box
use psi6_module
implicit none

real(dp), parameter :: sigma0 = 1._dp 
type(particledat), dimension(:), pointer :: particles
real(dp) :: cutoff
integer :: io_status
!type(cylformatter) :: cf
real(dp) :: p2
real(dp), dimension(3) :: direction
type(poly_box) :: simbox
type(factory) :: afactory
integer, parameter :: stdin = 5
!! The loop goes through all configurations fed to standard input of this 
!! program.
!cf = new_cylformatter()
cutoff = 2._dp * sigma0
do  
  !call readstate(cf, simbox, particles, io_status)
  call readstate(afactory, stdin, simbox, particles, io_status)
  if (io_status < 0) exit
  call globalnormal(simbox, pack(particles, particles%rod), cutoff, p2, direction)
  write(*,'(5('// fmt_char_dp() // ',1X))') psi6_bulk(simbox, pack(particles, particles%rod), direction), p2, direction
end do

end program
