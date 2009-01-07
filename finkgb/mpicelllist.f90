module mpicelllist
use mpi
use linkedlist
use particle, only : particledat

type pp 
  type(listtype), pointer :: list
end type pp

integer, save :: nprocs
integer, save :: thisproc,myix,myiy
integer, dimension(:,:), pointer,save :: procmatrix
real(dp),save :: cellside,domside
type(pp), dimension(:,:), pointer, save :: HOCs
integer, save :: ixlow,ixup,iylow,iyup
real(dp), save :: xlow,xup,ylow,yup


 
contains 

subroutine initCellList(parray,np,domw,cutoff,maxdr,xlow0,xup0,ylow0,yup0,&
                       &radius)
implicit none
type(particledat),dimension(:),intent(in) :: parray  
integer,intent(in) :: np, domw
real(dp), intent(in) :: cutoff, maxdr, xlow0, xup0, ylow0, yup0, radius
integer :: procsneeded
integer :: rc 
real(dp) :: mindomside
  xlow=xlow0
  ylow=ylow0
  xup=xup0
  xlow=xlow0
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,rc)
  call MPI_COMM_RANK(MPI_COMM_WORLD,thisproc,rc)
  mindomside=domw*(cutoff+2*maxdr)
  call decomposeCircular(radius,mindomside,procsneeded)
  if (procsneeded/=nprocs) stop 'initCellList: wrong number of ranks!'
  cellside=domside/domw
  ixlow=(myix-1)*domw; 
  ixup=myix*domw+1; 
  iylow=(myiy-1)*domw; 
  iyup=myiy*domw+1; 
  call newCellList(parray,np)
end subroutine initCellList


subroutine myParticles(myparray, mynp)
implicit none
type(particledat),dimension(:),pointer :: myparray, listparray
integer,intent(out) :: mynp
integer :: k, listnp, ix, iy, np0
np0=1000
myparray=>allocateParticleArray(np0)
k=1
do ix=lbound(HOCs,1)+1,ubound(HOCs,1)-1
do iy=lbound(HOCs,2)+1,ubound(HOCs,2)-1
  call toArray(HOCs(ix,iy)%list,listparray,listnp)
  if(listnp==0) cycle
  if((k+listnp)>np0) then
    myparray=>reallocateParticleArray(myparray,(np0+2*listnp))
    np0=np0+2*listnp
  end if
  myparray(k:k+listnp-1)=listparray(1:listnp)
  k=k+listnp
  !write(*,*) 'tehtiin taulukko solusta', ix, iy
end do
end do
mynp=k-1
end subroutine myParticles


subroutine newCellList(parray,np)
implicit none
type(particledat), dimension(:), intent(in) :: parray
integer, intent(in) :: np
integer :: astat,ix,iy,i
type(listtype), pointer :: head
real(dp) :: xmin, ymin
xmin=0
ymin=0 
  if(associated(HOCs)) then
    do ix=ixlow,ixup
    do iy=iylow,iyup
      call freeList(HOCs(ix,iy)%list)
    end do
    end do
  end if
  allocate(HOCs(ixlow:ixup,iylow:iyup),stat=astat)  
  if(astat/=0) stop 'newCellList: could not allocate HOCs table'
  do ix=lbound(HOCs,1),ubound(HOCs,1)
  do iy=lbound(HOCs,2),ubound(HOCs,2)
    nullify(HOCs(ix,iy)%list)
  end do
  end do
  do i=1,np
    ix=ceiling((parray(i)%x-xlow)/cellside)
    if(ix<ixlow .or. ix>ixup) cycle
    iy=ceiling((parray(i)%y-ylow)/cellside)
    if(iy<iylow .or. iy>iyup) cycle
    call prependToList(HOCs(ix,iy)%list,parray(i))
    if(ix>5) write(*,*) parray(i)
  end do
end subroutine newCellList



subroutine decomposeCircular(radius,mindomside,procsneeded)
use nrtype
implicit none
real(dp),intent(in) :: radius
real(dp),intent(in) :: mindomside
integer,intent(out) :: procsneeded
integer :: k,i,j,astat
integer :: nb1,domsperside
real(dp) :: b1,Rc
complex(dpc) :: z1,z2,z3,z4
  domsperside=floor((xup-xlow)/mindomside)
  nb1=domsperside
  b1=(xup-xlow)/nb1
  domside=b1
  Rc=radius
  k=0
  allocate(procmatrix(nb1,nb1), stat=astat)
  if(astat/=0) stop 'divide: procmatrix could not be allocated' 
  do i=1,nb1
    do j=1,nb1
      z1=cmplx((i-1)*b1-Rc,(j-1)*b1-Rc,dpc)
      z2=cmplx((i-1)*b1-Rc,j*b1-Rc,dpc)
      z3=cmplx(i*b1-Rc,(j-1)*b1-Rc,dpc)
      z4=cmplx(i*b1-Rc,j*b1-Rc,dpc)
      if (abs(z1)<Rc .or. abs(z2)<Rc .or. abs(z3)<Rc .or. abs(z4)<Rc) then
        procmatrix(i,j)=k
        if (thisproc==k) then
          myix=i
          myiy=j
        end if 
        k=k+1
      else 
        procmatrix(i,j)=-1
      end if
    end do
  end do
  procsneeded=k
end subroutine decomposeCircular


end module mpicelllist
