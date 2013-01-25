!> Reads configurations of Gay-Berne particles from a file. For each 
!! calculates and prints the tau1 translational ordering parameter. 
!!
!! Usage: program < configrationfile
!!
program tau1
  use nrtype
  use particle, only : particledat
  use tau1_module, only: tau1_routine
  use orientational_ordering, only: orientation_parameter
  use utils
  use class_factory
  implicit none
  type(particledat), dimension(:), pointer :: particles
  integer :: io_status
  real(sp) :: value
  real(sp) :: layer_distance
  real(dp), dimension(3) :: direction
  real(dp) :: P2
  type(factory) :: afactory
  integer, parameter :: stdin = 5
  type(poly_box) :: simbox
  do  
    call readstate(afactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else
      call orientation_parameter(particles, size(particles), P2, direction)
      call tau1_routine(particles, real(direction, sp), value, layer_distance)
      write(*, '(6(' // fmt_char_dp() // ',1X))') value, layer_distance, p2, direction
    end if
  end do
  if(associated(particles)) deallocate(particles)
end program 
