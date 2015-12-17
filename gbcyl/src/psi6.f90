!> Reads the given file configuration by configuration. For each configuration
!! calculates and prints the value of the psi6 with respect to the LC director,
!! the orientational ordering parameter and the director. 
!!
!! Usage: program < configuration_file
!!
program total_psi6
  use num_kind
  use particle, only : particledat
  use orientational_ordering
  use psi6_module, only : psi6_bulk
  use class_poly_box, only: poly_box
  use utils, only: fmt_char_dp
  use m_particle_factory, only: factory, factory_readstate
  implicit none
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  real(dp), dimension(3) :: direction
  real(dp) :: value
  integer, parameter :: stdin = 5
  type(poly_box) :: simbox
  type(factory) :: afactory
  do  
    call factory_readstate(afactory, stdin, simbox, particles, io_status)  
    if (io_status < 0) then
      exit
    else
       call orientation_parameter(pack(particles, particles%rod), &
            count(particles%rod), value, direction)
       write(*,'(5(' // fmt_char_dp() //',1X))') psi6_bulk(simbox, &
            pack(particles, particles%rod), direction), value, direction
    end if
  end do
end program


