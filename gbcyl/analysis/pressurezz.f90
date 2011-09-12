program pressurezz
use mc_engine
use class_factory
use particle
use class_poly_box
use class_parameterizer
use m_fileunit
use class_poly_nbrlist
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
real(dp) :: oldenergy, newenergy, dUinc, dUdec
real(dp) :: pzzinc, pzzdec, pzz
logical :: overlap
integer, parameter :: ndlz = 3
real(dp), dimension(ndlz) :: dLzs = (/0.001_dp, 0.003_dp, 0.01_dp/)
real(dp) :: dlz
integer :: i
character(len=3) :: idchar = '0'
character(len=1) :: direction = 'z'
type(poly_nbrlist) :: nbrlist
real(dp) :: newvolume, oldvolume, dVdec, dVinc
real(dp) :: temperature
read(*, *) idchar, direction
reader = new_parameterizer('parameters.in.' // trim(adjustl(idchar)))
call energy_init(reader)
call pnl_init(reader)
call getparameter(reader, 'temperature', temperature)
coordinateunit = fileunit_getfreeunit()
open(unit=coordinateunit, file='configurations.' //trim(adjustl(idchar)), action='READ', status='OLD')
write(*, '(6(A6,18X))') '#dV   ', 'P ', 'dV ', 'P', 'dV ', 'P '
do 
  call readstate(coordinatereader, coordinateunit, simbox, particles, ios)
  if (ios /= 0) then
    exit
  end if
  nparticles = size(particles)
  nbrlist = create_nbrlist(simbox, particles)
  !! Calculate initial volume and potential energy
  call potentialenergy(simbox, particles, nbrlist, oldenergy, overlap)
  oldvolume = volume(simbox)

  !! Start repeat for different dLz
  do i=1, ndlz
    dlz =dlzs(i)
    !! perform a increase in volume
    if(direction == 'z') then
      particles%z = particles%z*(getz(simbox) + dLz)/getz(simbox)
      call setz(simbox, getz(simbox) + dLz)
    else if(direction == 'y') then
      particles%y = particles%y*(gety(simbox) + dLz)/gety(simbox)
      call sety(simbox, gety(simbox) + dLz)
    else if(direction == 'x') then
      particles%x = particles%x*(getx(simbox) + dLz)/getx(simbox)
      call setx(simbox, getx(simbox) + dLz)
    end if
    !! Record new potential energy and volume
    call update(nbrlist, simbox, particles)
    call potentialenergy(simbox, particles, nbrlist, newenergy, overlap)
    if (overlap) then
      stop 'Overlap in trial volume change'
    end if
    newvolume = volume(simbox)
    !! record dV/V and change in energy dU
    dUinc = newenergy - oldenergy
    dVinc = newvolume - oldvolume
    pzzinc = temperature/dVinc*(real(nparticles,dp)*log(newvolume/oldvolume)-&
      dUinc/temperature)  

    !! Return old volume
    if(direction == 'z') then
      particles%z = particles%z*(getz(simbox) - dLz)/getz(simbox)
      call setz(simbox, getz(simbox) - dLz)
    else if(direction == 'y') then
      particles%y = particles%y*(gety(simbox) - dLz)/gety(simbox)
      call sety(simbox, gety(simbox) - dLz)
    else if(direction == 'x') then
      particles%x = particles%x*(getx(simbox) - dLz)/getx(simbox)
      call setx(simbox, getx(simbox) - dLz)
    end if

    !! perform a decrease in volume
    if (direction == 'z') then
      particles%z = particles%z*(getz(simbox) - dLz)/getz(simbox)
      call setz(simbox, getz(simbox) - dLz)
    else if(direction == 'y') then
      particles%y = particles%y*(gety(simbox) - dLz)/gety(simbox)
      call sety(simbox, gety(simbox) - dLz)
    else if(direction == 'x') then
      particles%x = particles%x*(getx(simbox) - dLz)/getx(simbox)
      call setx(simbox, getx(simbox) - dLz)
    end if

    !! Record new potential energy and volume
    call update(nbrlist, simbox, particles)
    call potentialenergy(simbox, particles, nbrlist, newenergy, overlap)
    if (overlap) then
      stop 'Overlap in trial volume change'
    end if
    newvolume = volume(simbox)
    !! record dV/V and change in energy dU
    dUdec = newenergy - oldenergy
    dVdec = newvolume - oldvolume
    pzzdec = temperature/dVdec*(real(nparticles,dp)*log(newvolume/oldvolume)-&
      dUdec/temperature)

    !! Calculate average pressure component from the increase and decrease of volume
    pzz = (pzzinc + pzzdec)/2._dp

    !! Write results to standard output
    write(*, '(2('// fmt_char_dp() //',1X))', ADVANCE='NO') dVinc, pzz


    !! Return old volume
    if(direction == 'z') then
      particles%z = particles%z*(getz(simbox) + dLz)/getz(simbox)
      call setz(simbox, getz(simbox) + dLz)
    else if(direction == 'y') then
      particles%y = particles%y*(gety(simbox) + dLz)/gety(simbox)
      call sety(simbox, gety(simbox) + dLz)
    else if(direction == 'x') then
      particles%x = particles%x*(getx(simbox) + dLz)/getx(simbox)
      call setx(simbox, getx(simbox) + dLz)
    end if

  end do !! End repeat for different dLz
  write(*, *) ''

  call delete(nbrlist)
end do

end program
