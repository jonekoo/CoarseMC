program doublezsize
use class_poly_box
use particle
use class_factory
implicit none

type(poly_box) :: simbox
type(particledat), allocatable :: particles(:)
type(factory) :: reader
type(factory) :: writer
type(particledat), allocatable :: newparticles(:)
integer :: ios
integer :: i,n
integer, parameter :: stdin=5
integer, parameter :: stdout=6
!! Read configuration in with cylformatter
call readstate(reader, stdin, simbox, particles, ios)
if (ios /= 0) then
  stop 'doublezsize: Error reading configuration.'
end if
n=size(particles)
allocate(newparticles(2*n))
newparticles(1:n)=particles(1:n)
newparticles(n+1:2*n)=particles(1:n)
do i=1,n
  newparticles(i)%z=newparticles(i)%z-0.5_dp*getz(simbox)    
  newparticles(i+n)%z=newparticles(i+n)%z+0.5_dp*getz(simbox)
end do
call setz(simbox, 2._dp*getz(simbox))
call writestate(writer, stdout, simbox, newparticles)
end program
