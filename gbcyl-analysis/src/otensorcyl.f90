program otensorcyl
use io
use nrtype
use particle, only : particledat
use order, only : eigens
use division, only : radialDivision
implicit none
character(len=50) :: buf
character(len=40) :: file
integer :: div
real(dp) :: Rc,Lz
type(particledat), dimension(:), pointer :: parray
real(dp), dimension(3,3) :: vectors
real(dp), dimension(3) :: diagonal
integer, dimension(:), pointer :: nps
type(particledat), dimension(:,:), pointer :: parrays
real(dp), dimension(:), pointer :: radia
integer :: np,i, astat
type(particledat), dimension(:), allocatable :: temp
  if (iargc()<1) stop 'Program needs arguments [file] [division]'
  call getarg(1,buf)
  read(unit=buf, fmt=*) file
  call getarg(2,buf) 
  read(unit=buf,fmt=*) div 
  call readstate(file,parray,Rc,Lz)
  np=size(parray)
  call radialDivision(Rc, div, parray, np, parrays, nps)
  allocate(temp(1:np), stat=astat)
  if(astat/=0) stop 'otensorcyl: allocation failed'
  do i=1,div
    temp(1:nps(i))=parrays(i,1:nps(i))    
    call eigens(temp,nps(i),diagonal,vectors)
    write(*,*) diagonal
  end do
  if(associated(parray)) deallocate(parray)
  if(associated(nps)) deallocate(nps)
  if(associated(parrays)) deallocate(parrays)
end program otensorcyl
