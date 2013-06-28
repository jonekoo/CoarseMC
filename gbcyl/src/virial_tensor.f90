program virial_tensor
!!
!! Calculates and prints the virial tensor for a group of Gay-Berne
!! and/or Lennard-Jones particles. The particles and the simulation box are 
!! read from the file configurations.id and parameters for the potential are 
!! read from file parameters.in.id. id is given to the program in standard 
!! input at the beginning of the program. 
!!
use class_factory !! for reading configurations.
use particle !! for particletype.
use class_poly_box !! for simulation box handling.
use class_parameterizer !! for reading parameters.
use m_fileunit !! to get a free output unit.
use class_pair_potential !! for forces
use nrtype !! for real kind=dp.
use utils !! for formatting output for real kind=dp.
implicit none
type(factory) :: coordinatereader
integer :: coordinateunit
integer :: nparticles
integer :: ios
type(particledat), allocatable :: particles(:)
type(poly_box) :: simbox
type(parameterizer) :: reader
character(len=3) :: idchar = '0'
real(dp) :: temperature
integer :: i, j, k, l
real(dp), dimension(3, 3) :: w
real(dp), dimension(3) :: rij, deriv
real(dp) :: cutoff
read(*, *) idchar
reader = new_parameterizer('inputparameters.'//trim(adjustl(idchar)))
!! Initialize the modules needed
call pp_init(reader)
call getparameter(reader, 'temperature', temperature)
call getparameter(reader, 'r_cutoff', cutoff)
coordinateunit = fileunit_getfreeunit()
open(unit=coordinateunit, file='configurations.' //trim(adjustl(idchar)), action='READ', status='OLD')
do
  w(:,:) = 0._dp 
  call readstate(coordinatereader, coordinateunit, simbox, particles, ios)
  if (ios /= 0) then
    exit
  end if
  nparticles = size(particles)
  do i=1, nparticles-1
    if (.not. particles(i)%rod) cycle
    do j=i+1, nparticles
      if (.not. particles(j)%rod) cycle
      !! Calculate the intermolecular position vector and
      rij = minimage(simbox, position(particles(i))-position(particles(j)))
       if (dot_product(rij, rij) > cutoff**2) cycle
      !! the gradient of the intermolecular potential (force)
      !forall(k=1:3) deriv(k) = d_potential(orientation(particles(i)), &
      !  orientation(particles(j)), rij, k)
      deriv = pair_force(particles(i), particles(j), rij)
      !! Calculate the virial tensor
      forall (k=1:3, l=1:3) w(k, l) = w(k, l) + rij(k) * deriv(l)
    end do
  end do
  w(:,:) = w(:,:)/volume(simbox)
  !forall(k=1:3) w(k, k) = w(k, k) + &
  !  real(nparticles,dp)*temperature/volume(simbox) 
  !  !! the above line may be sensitive to changes in epsilon0
  write(*, '(9('// fmt_char_dp() //',1X))') w(1, :), w(2, :), w(3, :)
end do


end program
