module particlearray
use nrtype 
use particle
implicit none

type(particledat), dimension(:), pointer :: thisparray

contains
  
  subroutine particleToXyz(N,particlearray,xs,ys,zs,uxs,uys,uzs,rods)
  implicit none
  type(particledat),dimension(:),pointer:: particlearray
  real(dp),dimension(:),pointer :: xs,ys,zs,uxs,uys,uzs
  logical,dimension(:),pointer :: rods
  integer,intent(in) :: N
  integer :: astat
    allocate(xs(N),ys(N),zs(N),uxs(N),uys(N),uzs(N),rods(N),stat=astat)
    if(astat/=0) stop 'particle_to_xyz: virhe muistia varattaessa'
    xs=particlearray%x
    ys=particlearray%y
    zs=particlearray%z
    uxs=particlearray%ux
    uys=particlearray%uy
    uzs=particlearray%uz
    rods=particlearray%rod
  end subroutine particleToXyz



  subroutine xyzToParticle(N,xs,ys,zs,uxs,uys,uzs,rods,particlearray)
  implicit none
  type(particledat),dimension(:),pointer :: particlearray 
  real(dp),dimension(:),intent(in) :: xs,ys,zs,uxs,uys,uzs
  logical,dimension(:),intent(in) :: rods
  integer,intent(in) :: N
  integer :: astat,i
    allocate(particlearray(N),stat=astat)
    if(astat/=0) stop 'xyz_to_particle: virhe muistia varattaessa'
    particlearray%x=xs
    particlearray%y=ys
    particlearray%z=zs
    particlearray%ux=uxs
    particlearray%uy=uys
    particlearray%uz=uzs
    particlearray%rod=rods
!    do i=1,N
!      particlearray(i)%index=i
!    end do
  end subroutine xyzToParticle


  function reallocateParticleArray(parray,np) result(newparray)
  implicit none
  type(particledat), dimension(:),pointer :: parray,newparray
  integer,intent(in) :: np
  integer :: astat,np0
    np0=size(parray)
    if(np0<np) then
      allocate(newparray(np),stat=astat)
      if(astat/=0) stop 'reallocateParticleArray: allocation failed'
      newparray(1:np0)=parray(1:np0)
      if(associated(parray)) deallocate(parray)
    else
      newparray=>parray(1:np)
    end if  
  end function reallocateParticleArray
 
  function allocateParticleArray(np) result(parray)
  implicit none
  integer,intent(in) :: np
  type(particledat),dimension(:),pointer :: parray
  integer :: astat
    allocate(parray(np), stat=astat)
    if (astat/=0) stop 'allocateParticleArray: allocation failed'
  end function allocateParticleArray


end module particlearray
