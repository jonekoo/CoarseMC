!! Calculates and prints the virial pressure tensor for a group of Gay-Berne
!! particles. The particles and the simulation box are read from the file 
!! configurations.id and parameters for the potential are read from file 
!! parameters.in.id. id is given to the program in standard input at the 
!! beginning of the program. 
program virialpressure
use class_factory
use particle
use class_poly_box
use class_parameterizer
use m_fileunit
use gayberne
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
real(dp) :: temperature
integer :: i, j, k, l
real(dp), dimension(3, 3) :: w
real(dp), dimension(3) :: rij, deriv 
read(*, *) idchar
reader = new_parameterizer('parameters.in.' // trim(adjustl(idchar)))
!! Initialize the modules needed
call gayberne_init(reader)
call getparameter(reader, 'temperature', temperature)
coordinateunit = fileunit_getfreeunit()
open(unit=coordinateunit, file='configurations.' //trim(adjustl(idchar)), action='READ', status='OLD')
do
  w(:,:) = 0._dp 
  !! Read particle coordinates
  call readstate(coordinatereader, coordinateunit, simbox, particles, ios)
  if (ios /= 0) then
    exit
  end if
  nparticles = size(particles)
  do i=1, nparticles-1
    do j=i+1, nparticles
      !! Calculate the intermolecular position vector and
      rij = minimage(simbox, position(particles(i))-position(particles(j)))
      !! the gradient of the intermolecular potential (force)
      forall(k=1:3) deriv(k) = d_potential(orientation(particles(i)), &
        orientation(particles(j)), rij, k)
      !! Calculate the virial tensor
      forall(k=1:3, l=1:3) w(k, l) = w(k, l) + rij(k)*deriv(l)
    end do
  end do
  w(:,:) = -w(:,:)/volume(simbox)
  forall(k=1:3) w(k, k) = w(k, k) + &
    real(nparticles,dp)*temperature/volume(simbox) 
  !! the above line may be sensitive to changes in epsilon0
  write(*, *) w(1, :)
  write(*, *) w(2, :)
  write(*, *) w(3, :)
  write(*, *) ''
end do


end program
