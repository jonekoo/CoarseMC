program printforces
use num_kind
use class_pair_potential
use utils, only : fmt_char_dp
use m_fileunit
use particlewall
use class_poly_box
use particle, only: particledat, setposition, setorientation
use m_particle_factory, only: factory, factory_readstate
use class_parameterizer
implicit none

real(dp), dimension(3), parameter :: ex = (/1._dp, 0._dp, 0._dp/)
real(dp), dimension(3), parameter :: ey = (/0._dp, 1._dp, 0._dp/)
real(dp), dimension(3), parameter :: ez = (/0._dp, 0._dp, 1._dp/)
character(4) :: idchar
logical :: is_wall_on

type(factory) :: coordinatereader
integer :: coordinateunit, ios
character(len=200) :: statefile
type(parameterizer) :: reader

type(particledat), allocatable :: particles(:)
type(poly_box) :: simbox

type(particledat) :: testparticles(4)
character(len=3) :: particle_descriptions(4)

integer :: i, j, k
integer, parameter :: n = 1000
type(particledat) :: temp
real(dp), parameter :: step = 0.01_dp
real(dp) :: r, rij(3)
real(dp) :: rmax
real(dp) :: f(3)
real(dp) :: energy
integer :: err

class(pair_interaction), allocatable :: pair_ia

read(*, *) idchar

!! Read inputparameters
reader = new_parameterizer('inputparameters.'//trim(adjustl(idchar)), &
     logfile = "printpotentials_log."//trim(adjustl(idchar)))
allocate(pair_ia, source=conditional_pair_interaction(reader))

call getparameter(reader, 'is_wall_on', is_wall_on)
if (is_wall_on) then
  call particlewall_init(reader)
end if

!! Read inputconfiguration
coordinateunit = fileunit_getfreeunit()
statefile = 'inputconfiguration.'//trim(adjustl(idchar))
!! In any case the output should be appended to configurations.(id) and restartfile
!! should be used only for reading once and then overwriting the file.
open(file=statefile, unit=coordinateunit, action='READ', status='OLD', iostat=ios)
call factory_readstate(coordinatereader, coordinateunit, simbox, particles, ios)
if (0/=ios) then 
  write(*, *) 'Error ', ios,' reading ', statefile, ' Stopping.' 
  stop
end if
close(coordinateunit)

!! :TODO: Create only the kinds of particle pairs that exist in the inputconfiguration
!! Look for GB particles:

j=1
do i = 1, size(particles)
  if (particles(i)%rod) then
    testparticles(j) = particles(i)
    call setorientation(testparticles(j), ex)
    particle_descriptions(j) = 'GBx'
    j = j + 1
    testparticles(j) = particles(i)
    call setorientation(testparticles(j), ey)
    particle_descriptions(j) = 'GBy'
    j = j + 1
    testparticles(j) = particles(i)
    call setorientation(testparticles(j), ez)
    particle_descriptions(j) = 'GBz'
    j = j + 1
    exit
  end if
end do

!! Look for LJ particles:
do i = 1, size(particles)
  if (.not. particles(i)%rod) then
    testparticles(j) = particles(i)
    particle_descriptions(j) = 'LJ'
    j = j + 1
    exit
  end if
end do

rmax = simbox%lx / 2._dp
do i = 1, size(testparticles)
  call setposition(testparticles(i), (/rmax, 0._dp, 0._dp/))
end do


write(*, '(A23,1X)', advance='no') 'r'
do j = 1, size(particle_descriptions)
   temp = testparticles(j)
   temp%x = testparticles(j)%x - r  
   if (is_wall_on) then
      write(*, '(A23,1X)', advance='no') trim(particle_descriptions(j))// &
           "-Wall"
   end if
   do k = 1, j
      write(*, '(A23,1X)', advance='no') trim(particle_descriptions(j)) // &
           "-" // particle_descriptions(k)
   end do
end do
write(*, '(A)') ''


r = 0._dp
!! Calculate all possible pair potentials:
do i = 1, n
  r = r + step
  if (r > rmax) exit
  write(*, fmt='(G23.15E3,1X)', advance='no') r
  do j = 1, size(testparticles)
     temp = testparticles(j)
     temp%x = testparticles(j)%x - r  
     if (is_wall_on) then
        call particlewall_potential(temp, simbox, energy, err)
        if (err /= 0) then 
           write(*, '(A23,1X)', advance='no') 'NaN'
        else
           rij = position(temp)
           f = particlewall_force(temp, simbox) 
           !! The force has been defined as f=-dU/dr=-dU/dr_w * dr_w/dr
           write(*, '(3(G23.15E3,1X))', advance='no') dot_product(rij, f) / &
                sqrt(dot_product(rij, rij))
        end if
     end if
     do k = 1, j
        call pair_ia%pair_potential(temp, testparticles(k), &
             position(testparticles(k))-position(temp), energy, err)
        if (err /= 0) then 
           write(*, '(A23,1X)', advance='no') 'NaN'
        else
           rij = position(temp) - position(testparticles(k))
           f = pair_ia%pair_force(testparticles(k), temp, rij)
           ! We want to print out the force acting on temp by
           ! testparticles(k), when testparticles(k) is at origin.
           write(*, '(G23.15E3,1X)', advance='no') dot_product(rij, f) / &
                sqrt(dot_product(rij, rij))
        end if
     end do
  end do
  write(*, '(A)') ''
end do
end program
