program onwallpressure
use class_factory
use particle
use class_poly_box
use class_parameterizer
use m_fileunit
use energy
use nrtype
use utils
implicit none
type(factory) :: coordinatereader
integer :: coordinateunit
integer :: nparticles
integer :: ios
type(particledat), dimension(:), pointer :: particles
type(poly_box) :: simbox
type(parameterizer) :: reader
character(len=3) :: idchar = '0'
real(dp), parameter :: dr=1.e-6_dp
real(dp) :: f01, f12, f, totv0, totv1, totv2
logical :: overlap
read(*, *) idchar
reader = new_parameterizer('inputparameters.'//trim(adjustl(idchar)))
!! Initialize the modules needed
call energy_init(reader)
coordinateunit = fileunit_getfreeunit()
open(unit=coordinateunit, file='configurations.' //trim(adjustl(idchar)), action='READ', status='OLD')
do
  call readstate(coordinatereader, coordinateunit, simbox, particles, ios)
  if (ios /= 0) then
    exit
  end if
  nparticles = size(particles)
  call totwallprtclV(simbox, particles, totv1, overlap)
  call setx(simbox, getx(simbox)+dr)
  call totwallprtclV(simbox, particles, totv2, overlap)
  call setx(simbox, getx(simbox)-2._dp*dr)
  call totwallprtclV(simbox, particles, totv0, overlap)
  f01=(totv1-totv0)/dr
  f12=(totv2-totv1)/dr
  f=0.5_dp*(f01+f12)
  write(*, *) f, -f/(getz(simbox)*4._dp*atan(1._dp)*getx(simbox))
end do



end program
