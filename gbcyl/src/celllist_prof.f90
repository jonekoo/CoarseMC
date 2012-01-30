program celllist_prof
use particle
use class_poly_box
use class_factory
use cell_energy
use cell, only: list
use class_simplelist
use nrtype
implicit none

type(particledat), dimension(:), pointer :: particles
type(poly_box) :: simbox
type(factory) :: configuration_reader
integer, parameter :: stdin=5
real(dp), parameter :: minlength=2.8_dp
real(dp) :: t0,t1,t2
type(list) :: cl
type(simplelist) :: sl
integer :: i
integer :: ios
logical, dimension(:), allocatable :: slmask, clmask
!read a particlefile
call readstate(configuration_reader, stdin, simbox, particles, ios)
allocate(slmask(size(particles)), clmask(size(particles)))
!create cell lists
call simplelist_init(minlength, iseven=.false.,updatethreshold=0.5_dp)
sl=new_simplelist(simbox,particles)
call cell_energy_init(minlength, iseven=.false.)
cl=new_celllist(simbox,particles)
!measure creation times
!create nbrmasks for all particles
call cpu_time(t0)
do i=1,size(particles)
  slmask=nbrmask(sl,simbox,particles,i)
end do
call cpu_time(t1)
do i=1,size(particles)
  clmask=nbrmask(cl,simbox,particles,i)
end do
do i=1,size(particles)
  clmask=nbrmask(cl,simbox,particles,i)
  slmask=nbrmask(sl,simbox,particles,i)
  write(*, *) all(slmask.eqv.clmask)
end do
call cpu_time(t2)
write(*,*) all(slmask .eqv. clmask)
write(*,*) t1-t0, t2-t1
!measure nbrmask creation times
!compare nbrmask creation times
end program
