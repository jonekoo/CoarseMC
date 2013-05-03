!> Reads the given file configuration by configuration. For each configuration
!! calculates and prints the value of the psi6 with respect to the LC director,
!! the orientational ordering parameter and the director. 
!!
!! Usage: program < configuration_file
!!
program psi6
  use nrtype, only: dp
  use particle, only : particledat
  use orientational_ordering
  use psi6_module, only : psi6_bulk
  use class_poly_box
  use utils
  use class_factory
  implicit none
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  real(dp), dimension(3) :: direction
  real(dp) :: value
  integer, parameter :: stdin = 5
  type(poly_box) :: simbox
  type(factory) :: afactory
  do  
    call readstate(afactory, stdin, simbox, particles, io_status)  
    if (io_status < 0) then
      exit
    else
      call orientation_parameter(pack(particles, particles%rod), count(particles%rod), value, direction)
      write(*,'(5(' // fmt_char_dp() //',1X))') psi6_bulk(simbox, pack(particles, particles%rod), direction), value, direction
    end if
  end do
end program


