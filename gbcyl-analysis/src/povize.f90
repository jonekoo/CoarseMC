program povize
use io, only : readstate, povout
use nrtype
use particle, only : particledat
implicit none
character(len=50) :: buf
character(len=40) :: file
integer :: div
real(dp) :: Rc,Lz
type(particledat), dimension(:),pointer :: parray
real(dp), dimension(3,3) :: vectors
real(dp), dimension(3) :: diagonal

  if (iargc()/= 1) stop 'Program needs argument [particlefile]'
  call getarg(1,buf)
  read(unit=buf, fmt=*) file
  call readstate(file,parray,Rc,Lz) 
  call povout(parray,Rc,Lz)
  if(associated(parray)) deallocate(parray)
end program povize
