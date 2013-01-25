SUBROUTINE tred2(a,d,e,novectors) 
USE nrtype; USE nrutil, ONLY : assert_eq, outerprod 
IMPLICIT NONE 
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a 
REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e 
LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors 
! Householder reduction of a real, symmetric, N 
! × N matrix a. On output, a is replaced 
! by the orthogonal matrix Q eﬀecting the transformation. d returns the diagonal elements 
! of the tridiagonal matrix, and e the oﬀ-diagonal elements, with e(1)=0. If the optional 
! argument no vectors is present, only eigenvalues are to be found subsequently, in which 
! case a contains no useful information on output. 
INTEGER(I4B) :: i, j, l, n 
REAL(SP) :: f, g, h, hh, scale 
REAL(SP), DIMENSION(size(a,1)) :: gg 
LOGICAL(LGT), SAVE :: yesvec = .true. 
n=assert_eq(size(a,1),size(a,2),size(d),size(e),'tred2') 
if(present(novectors)) yesvec = .not. novectors 
do i=n,2,-1 
  l=i-1 
  h=0.0 
  if(l> 1)then 
    scale=sum(abs(a(i,1:l))) 
    if(scale==0.0)then ! Skip transformation. 
      e(i)=a(i,l) 
    else 
      a(i,1:l)=a(i,1:l)/scale ! Use scaled a’s for transformation. 
      h=sum(a(i,1:l)**2) ! Form σ in h. 
      f=a(i,l) 
      g=-sign(sqrt(h),f) 
      e(i)=scale*g 
      h=h-f*g !! Now h is equation (11.2.4). 
      a(i,l)=f-g !! Store u in the ith row of a. 
      if(yesvec)a(1:l,i)=a(i,1:l)/h !! Store u/H in ith column of a. 
      do j=1,l !! Store elements of p in temporarily unused elements of e. 
        e(j)=(dot_product(a(j,1:j),a(i,1:j))& 
        +dot_product(a(j+1:l,j),a(i,j+1:l)))/h 
      end do 
      f=dot_product(e(1:l),a(i,1:l)) 
      hh=f/(h+h) !! Form K, equation (11.2.11). 
      e(1:l)=e(1:l)-hh*a(i,1:l) 
      !! Form q and store in e overwriting p. 
      do j=1,l !! Reduce a, equation (11.2.13). 
        a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j) 
      end do 
    endif 
  else 
    e(i)=a(i,l) 
  endif 
  d(i)=h 
enddo 
if(yesvec) d(1)=0.0 
e(1)=0.0 
do i=1,n !! Begin accumulation of transformation matrices. 
  if(yesvec)then 
    l=i-1 
    if(d(i)/=0.0)then 
      !! This block skipped when i=1. Use u and u/H stored in a to form P·Q. 
      gg(1:l)=matmul(a(i,1:l),a(1:l,1:l)) 
      a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l)) 
    endif 
    d(i)=a(i,i) 
    a(i,i)=1.0 !! Reset row and column of a to identity matrix for next iteration. 
    a(i,1:l)=0.0 
    a(1:l,i)=0.0 
  else 
    d(i)=a(i,i) 
  endif 
enddo 
END SUBROUTINE tred2 

