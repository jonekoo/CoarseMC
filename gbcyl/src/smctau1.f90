!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Calculates and prints the translational tau1 parameter for a Sm-C or Sm-A 
!! phase. Other information printed is the layer distance, the orientational
!! ordering parameter of the local layer normalss and the reference direction,
!! the global layer normal used in the calculation of the tau1 parameter.
!!
!! Usage: smctau1 < configurationfile
!!
program smctau1

!! These are needed for the types read from configuration files
use particle, only: particledat
use class_poly_box, only: poly_box

!! For printing
use utils, only: fmt_char_dp

!! For type definitions
use num_kind

!! For reading configurations from disk
use m_particle_factory, only: factory, factory_readstate

!! To determine the normal of the layers
use layernormal

!! To calculate the smectic order parameter tau1
use tau1_module
implicit none

real(dp) :: sigma0 = 1._dp 
!real(dp), dimension(:, :), allocatable :: localnormals
type(particledat), allocatable :: particles(:)
!integer :: nparticles
real(dp) :: cutoff
!real(dp) :: rcyl, lz
integer :: io_status
!!type(cylformatter) :: cf
real(dp) :: p2
real(dp), dimension(3) :: direction
real(sp) :: tau1_value
real(sp) :: layerdistance
type(poly_box) :: simbox
type(factory) :: afactory
integer, parameter :: stdin = 5

!! The loop goes through all configurations fed to standard input of this 
!! program.
!!cf = new_cylformatter()
do  
  call factory_readstate(afactory, stdin, simbox, particles, io_status)
  if (io_status < 0) exit
  !if (allocated(localnormals) .and. size(localnormals)/3 /= nparticles) then
  !  deallocate(localnormals)
  !end if
  !if (.not. allocated(localnormals)) allocate(localnormals(nparticles, 3))
  !simbox = new_cylinder(rcyl, lz)
  cutoff = 2._dp * sigma0
  call globalnormal(simbox, pack(particles, particles%rod), cutoff, p2, direction)
  call tau1_routine(pack(particles, particles%rod), real(direction, sp), tau1_value, layerdistance)
  write(*, '(6('//fmt_char_dp()//',1X))') tau1_value, layerdistance, p2, direction
end do

end program
