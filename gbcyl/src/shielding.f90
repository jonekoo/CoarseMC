program shielding
  use nrtype, only: dp
  use particle, only : particledat
  use orientational_ordering
  use class_poly_box
  use utils
  use class_factory
  use m_shielding
  use class_parameterizer
  use m_fileunit
  use gblj
  implicit none
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  type(poly_box) :: simbox
  type(factory) :: afactory
  real(dp) :: tensor(3, 3)
  real(dp) :: values(3), vectors(3, 3)
  type(parameterizer) :: reader
  character(4) :: idchar
  integer :: coordinateunit
  read(*, *) idchar
  reader = new_parameterizer('inputparameters.'//trim(adjustl(idchar)), logfile = "virial_tensor_log."//trim(adjustl(idchar)))
  !! Initialize the modules needed
  call gblj_init(reader)
  coordinateunit = fileunit_getfreeunit()
  open(unit=coordinateunit, file='configurations.' //trim(adjustl(idchar)), action='READ', status='OLD')
  do  
    call readstate(afactory, coordinateunit, simbox, particles, io_status)  
    if (io_status < 0) then
      exit
    else
      call eigens(pack(particles, particles%rod), count(particles%rod), values, vectors)
      call cycle_largest_to_3rd(values, vectors)
      call avg_shielding(simbox, particles, axes = vectors, tensor = tensor)
      write(*, '( 9('// fmt_char_dp() //',1X))') tensor
    end if
  end do
end program


