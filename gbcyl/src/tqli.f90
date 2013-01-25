SUBROUTINE tqli(d,e,z) 
USE nrtype; USE nrutil,ONLY: assert_eq, nrerror 
USE nr, ONLY: pythag 
IMPLICIT NONE 
REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e 
REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT)::z 
!! QL algorithm with implicit shifts, to determine the eigenvalues and 
!! eigenvectors of a real, symmetric, tridiagonal matrix, or of a real, 
!! symmetric matrix previously reduced by tred2 §11.2. d is a vector of length
!! N. On input, its elements are the diagonal elements of the 
!! tridiagonal matrix. On output, it returns the eigenvalues. The vector e 
!! inputs the subdiagonal elements of the tridiagonal matrix, with e(1) 
!! arbitrary. On output e is destroyed. When ﬁnding only the eigenvalues,
!! the optional argument z is omitted. If the eigenvectors of a tridiagonal
!! matrix are desired, the N×N matrix z is input as the identity matrix. If 
!! the eigenvectors of a matrix that has been reduced by tred2 are required, 
!! then z is input as the matrix output by tred2. In either case, the kth 
!! column of z returns the normalized eigenvector corresponding to d(k). 
INTEGER(I4B)::i,iter,l,m,n,ndum 
REAL(SP)::b,c,dd,f,g,p,r,s 
REAL(SP),DIMENSION(size(e))::ff 
n=assert_eq(size(d),size(e),'tqli:n') 
if(present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli:ndum') 
e(:)=eoshift(e(:),1) !! Convenient to renumber the elements of e. 
do l=1,n 
  iter=0 
  iterate:do 
    do m=l,n-1 !! Look for a single small subdiagonal element to split the matrix. 
      dd=abs(d(m))+abs(d(m+1)) 
      if(abs(e(m))+dd==dd)exit 
    end do 
    if(m==l) exit iterate 
    if(iter==30) call nrerror('too many iterations in tqli') 
    iter=iter+1 
    g=(d(l+1)-d(l))/(2.0_sp*e(l)) !! Form shift. 
    r=pythag(g,1.0_sp) 
    g=d(m)-d(l)+e(l)/(g+sign(r,g)) !! This is d_m−k_s. 
    s=1.0 
    c=1.0 
    p=0.0 
    do i=m-1,l,-1 !! A plane rotation as in the original QL, followed by Givens
                  !! rotations to restore tridiagonal form. 
      f=s*e(i) 
      b=c*e(i) 
      r=pythag(f,g) 
      e(i+1)=r 
      if(r==0.0)then !! Recover from underﬂow. 
        d(i+1)=d(i+1)-p 
        e(m)=0.0 
        cycle iterate 
      endif 
      s=f/r 
      c=g/r 
      g=d(i+1)-p 
      r=(d(i)-g)*s+2.0_sp*c*b 
      p=s*r 
      d(i+1)=g+p 
      g=c*r-b 
      if(present(z))then !! Form eigenvectors. 
        ff(1:n)=z(1:n,i+1) 
        z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n) 
        z(1:n,i)=c*z(1:n,i)-s*ff(1:n) 
      endif 
    enddo 
    d(l)=d(l)-p 
    e(l)=g 
    e(m)=0.0 
  end do iterate 
enddo 
END SUBROUTINE tqli 
