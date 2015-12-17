!> Reads configurations of Gay-Berne particles from a file. For each 
!! calculates and prints the tau1 translational ordering parameter. 
!!
!! Usage: program < configrationfile
!!
program tau1_program
  use num_kind
  use particle, only : particledat
  use tau1_module, only: tau1_routine
  use orientational_ordering, only: orientation_parameter
  use utils
  use m_particle_factory
  implicit none
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  real(sp) :: value
  real(sp) :: layer_distance
  real(dp), dimension(3) :: direction
  real(dp) :: P2
  type(factory) :: afactory
  integer, parameter :: stdin = 5
  type(poly_box) :: simbox
  do  
    call factory_readstate(afactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else
       call orientation_parameter(pack(particles, particles%rod), &
            count(particles%rod), P2, direction)
       call tau1_routine(pack(particles, particles%rod), real(direction, sp), &
            value, layer_distance)
       write(*, '(6(' // fmt_char_dp() // ',1X))') value, layer_distance, p2, &
            direction
    end if
  end do
end program 
