module cutgeometry
use particle
use nrtype
implicit none

contains

subroutine readinatoms(N, Lx, Ly, Lz, particles)
  use nrtype, only: dp
  use particle, only: particledat
  implicit none

  type(particledat), dimension(:), intent(out) :: particles
  real(dp),intent(out) :: Lx,Ly,Lz
  integer, intent(out) :: N
  integer :: ios,i
  integer,dimension(:),allocatable :: help

  open(20,file='atoms.in',status='old',form='formatted',iostat=ios)
  if(ios/=0)then
     write(*,*)'Virhe tiedoston atoms.in avaamisessa. Ohjelma keskeytyy..'
     stop;
  end if
  read(20,*)N,Lx,Ly,Lz
  allocate(help(N),stat=ios)
  if(ios/=0)then
     write(*,*)'Virhe muistinvaraamisessa.Ohjelma keskeytyy..'
     stop;
  end if
  read(20,*)particles(1:N)%x,particles(1:N)%y,particles(1:N)%z,particles(1:N)%ux,particles(1:N)%uy,particles(1:N)%uz,help(1:N)
  close(20);

  do i=1,N
     if(help(i)==1)particles(i)%rod=.true.
     if(help(i)==0)particles(i)%rod=.false.
  end do
end subroutine readinatoms



!! Forms a new array of particles by discarding particles outside of 
!! @p radius. 
!!
!! param particles the array of particles.
!! param n_particles number of particles in @p particles.
!! param radius the radius used to discard or accept particles.
!! param new_particles new array of particles.
!! param n_new_particles number of particles in @p particles.
!!
subroutine cut_with_radius(particles, n_particles, radius, cut, & 
new_particles, n_new_particles)
  type(particledat), dimension(:), pointer :: particles
  integer, intent(in) :: n_particles
  real(dp), intent(in) :: radius
  real(dp), intent(in) :: cut
  type(particledat), dimension(:), pointer :: new_particles
  integer, intent(out) :: n_new_particles
  type(particledat), dimension(:), allocatable :: tempptarray
  integer :: i, j = 0, astat
  real(dp) :: radius_i
  allocate(tempptarray(n_particles), stat = astat)
  if(astat/=0) then 
    write (*,*) 'initptarray: Virhe varattaessa muistia tilapäiselle'
    write (*,*) 'hiukkastaulukolle. Ohjelma keskeytyy.'
    stop;
  end if
  do i = 1, n_particles
    radius_i = particles(i)%x*particles(i)%x + particles(i)%y*particles(i)%y
    radius_i = sqrt(radius_i)
    if (radius-radius_i > cut ) then
      j = j + 1
      tempptarray(j) = particles(i)
    end if
  end do
  n_new_particles = j
  allocate(new_particles(n_new_particles), stat = astat)
  if(astat /= 0) then
    write(*,*) 'initptarray: Virhe varattaessa muistia hiukkastaulukolle.'
    write(*,*) 'Ohjelma keskeytyy.'
    stop;
  end if
  new_particles = tempptarray(1:n_new_particles)
  if (allocated(tempptarray)) deallocate(tempptarray)
end subroutine 

!! Forms a new array of particles by discarding particles that don't fill condition.
!!
!! @p particles the array of particles.
!! @p n_particles number of particles in @p particles.
!! @p condition the condition function to accept the particles
!! @p new_particles new array of particles.
!! @p n_new_particles number of particles in @p particles.
!!
subroutine cut_with_condition(particles, n_particles, condition, new_particles, &
n_new_particles)
  implicit none
  interface
    function condition(a_particle) result(istrue)
      use particle
      type(particledat), intent(in) :: a_particle
      logical :: istrue
    end function 
  end interface
  type(particledat), dimension(:), pointer :: particles
  integer, intent(in) :: n_particles
  type(particledat), dimension(:), pointer :: new_particles
  integer, intent(out) :: n_new_particles
  type(particledat), dimension(:), allocatable :: tempptarray
  integer :: i, j = 0, astat
  allocate(tempptarray(n_particles), stat = astat)
  if(astat/=0) then 
    write (*,*) 'initptarray: Virhe varattaessa muistia tilapäiselle'
    write (*,*) 'hiukkastaulukolle. Ohjelma keskeytyy.'
    stop;
  end if
  do i = 1, n_particles
    if (condition(particles(i))) then
      j = j + 1
      tempptarray(j) = particles(i)
    end if
  end do
  n_new_particles = j
  allocate(new_particles(n_new_particles), stat = astat)
  if(astat /= 0) then
    write(*,*) 'initptarray: Virhe varattaessa muistia hiukkastaulukolle.'
    write(*,*) 'Ohjelma keskeytyy.'
    stop;
  end if
  new_particles = tempptarray(1:n_new_particles)
  if (allocated(tempptarray)) deallocate(tempptarray)
end subroutine

end module

program tightpack
  use nrtype
  use particle
  use class_cylformatter
  use cutgeometry
  implicit none 

  integer :: Nalloc = 104544, astat
  type(particledat), dimension(:), pointer :: particles
  integer :: n_particles
  type(particledat), dimension(:), pointer :: new_particles
  integer :: n_new_particles
  real(dp) :: Lx, Ly, Lz
  real(dp) :: radius
  real(dp) :: offset
  character(len = *), parameter :: outfilename = 'cutcylinder.out'
  type(cylformatter) :: cf

  write(*, *) 'Radius: ' 
  read(*, *) radius 
  write(*, *) 'Offset (positive) from the radius: ' 
  read(*, *) offset 
  allocate(particles(Nalloc), stat = astat)
  if(astat /= 0) stop 'tightpack.f90: Virhe varattaessa muistia'
  call readinatoms(n_particles, Lx, Ly, Lz, particles)
  call cut_with_radius(particles, n_particles, radius, offset,&
                     & new_particles, n_new_particles)
  deallocate(particles)  
  cf = new_cylformatter(outfilename)
  call writestate(cf, new_particles, n_new_particles, radius, Lz)
  deallocate(new_particles)

end program tightpack


