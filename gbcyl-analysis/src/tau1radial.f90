program tau1radial
use io
use nrtype
use particle, only : particledat
use order, only : tau1
implicit none
character(len=50) :: buf
character(len=40) :: file
real(dp) :: Rc,Lz,d
type(particledat), dimension(:), pointer :: parray
real(dp), dimension(:), pointer :: rs
real(dp), dimension(3) :: dirvec
integer :: np,i, astat
  if (iargc()<1) stop 'Program needs arguments [file] [distance between layers]'
  call getarg(1,buf)
  read(unit=buf, fmt=*) file
  call getarg(2,buf)
  read(unit=buf, fmt=*) d
  call readstate(file,parray,Rc,Lz)
  np=size(parray)
  allocate(rs(1:np), stat=astat)
  if(astat/=0) stop 'tau1radial: memory allocation failed'
  do i=1,np
    rs(i)=sqrt( (parray(i)%x)**2 + (parray(i)%y)**2) 
  end do
  write(*,*) tau1(rs, np, d) 
  if(associated(parray)) deallocate(parray)
  if(associated(rs)) deallocate(rs)
end program tau1radial
