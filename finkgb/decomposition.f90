module decomposition
implicit none

!! Cells per one domainside
integer, save :: cellsperside 

!! This process's indices in the processmatrix
integer, save :: myix,myiy 

!! Processmatrix dimensions
integer, save :: matdimx, matdimy 


contains  

!! Uses the rule cutter to make a matrix of processes.
!! Needs arguments
!! b1=cell sidelength
!! nb1=number of cells in one direction
!! cutter=user supplied function
!! Returns
!! procmatrix=matrix of processes 
!! nproc=number of processes
subroutine divide(b1,nb1,Rc,procmatrix,nproc)
use nrtype
implicit none
integer,intent(in) :: nb1
integer,dimension(:,:),pointer :: procmatrix
integer,intent(out) :: nproc
real(dp),intent(in) :: b1,Rc
integer :: k,i,j
complex(dpc) :: z1,z2,z3,z4
  matdimx=nb1
  matdimy=nb1
  k=0
  do i=1,nb1
    do j=1,nb1
      z1=cmplx((i-1)*b1-Rc,(j-1)*b1-Rc,dpc)
      z2=cmplx((i-1)*b1-Rc,j*b1-Rc,dpc)
      z3=cmplx(i*b1-Rc,(j-1)*b1-Rc,dpc)
      z4=cmplx(i*b1-Rc,j*b1-Rc,dpc)
      if (abs(z1)<Rc .or. abs(z2)<Rc .or. abs(z3)<Rc .or. abs(z4)<Rc) then
         procmatrix(i,j)=k
         k=k+1
      else 
        procmatrix(i,j)=-1
      end if
    end do
  end do
  nproc=k
end subroutine divide



function Rcutter(x,y,z,params)
  use nrtype
  implicit none
  logical :: Rcutter
  real(dp) :: x,y,z
  real(dp),dimension(:) :: params
  if(x**2+y**2<params(1)**2) then 
    Rcutter = .true.   
  else 
    Rcutter = .false.
  end if
end function Rcutter





end module decomposition
