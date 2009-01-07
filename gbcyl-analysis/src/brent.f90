FUNCTION brent(ax,bx,cx,func,tol,xmin) 
USE nrtype; USE nrutil, ONLY : nrerror 
IMPLICIT NONE 
REAL(SP), INTENT(IN) :: ax,bx,cx,tol 
REAL(SP), INTENT(OUT) :: xmin 
REAL(SP) ::brent 
INTERFACE 
  FUNCTION func(x) 
  USE nrtype 
  IMPLICIT NONE 
  REAL(SP), INTENT(IN) :: x 
  REAL(SP) :: func 
  END FUNCTION func 
END INTERFACE 
INTEGER(I4B), PARAMETER :: ITMAX=100 
REAL(SP), PARAMETER :: CGOLD=0.3819660_sp,ZEPS=1.0e-3_sp*epsilon(ax) 
  !! Given a function func, and given a bracketing triplet of abscissas ax, bx,
  !! cx (such that bx is between ax and cx, and func(bx) is less than both 
  !! func(ax) and func(cx)), this routine isolates the minimum to a fractional
  !! precision of about tol using Brent’s method. 
  !! The abscissa of the minimum is returned as xmin, and the minimum function 
  !! value is returned as brent, the returned function value. 
  !! Parameters: Maximum allowed number of iterations; golden ratio; and a 
  !! small number that protects against trying to achieve fractional accuracy 
  !! for a minimum that happens to be exactly zero. 
INTEGER(I4B) :: iter 
REAL(SP) ::a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm 
!! a and b must be in ascending order, though 
!! the input abscissas need not be. 
a=min(ax,cx) 
b=max(ax,cx) 
v=bx !! Initializations... 
w=v 
x=v 
e=0.0 !! This will be the distance moved on the step before last. 
fx=func(x) 
fv=fx 
fw=fx 
do iter=1,ITMAX !! Main program loop. 
  xm=0.5_sp*(a+b) 
  tol1=tol*abs(x)+ZEPS 
  tol2=2.0_sp*tol1 
  if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) then !! Test for done here. 
    xmin=x !! Arrive here ready to exit with best values. 
    brent=fx 
    RETURN 
  end if 
  if (abs(e) > tol1) then !! Construct a trial parabolic ﬁt. 
    r=(x-w)*(fx-fv) 
    q=(x-v)*(fx-fw) 
    p=(x-v)*q-(x-w)*r 
    q=2.0_sp*(q-r) 
    if (q> 0.0) p=-p 
    q=abs(q) 
    etemp=e 
    e=d 
    if (abs(p) >= abs(0.5_sp*q*etemp) .or. & 
      p<= q*(a-x) .or. p>= q*(b-x)) then 
      !! The above conditions determine the acceptability of the parabolic ﬁt.
      !! Here it is not o.k.,so we take the golden section step into the 
      !! larger of the two segments. 
      e=merge(a-x,b-x, x >= xm ) 
      d=CGOLD*e 
    else !! Take the parabolic step. 
      d=p/q 
      u=x+d 
      if(u-a <tol2 .or. b-u < tol2) d=sign(tol1,xm-x) 
    end if 
  else !! Take the golden section step into the larger of the two segments. 
    e=merge(a-x,b-x, x>= xm ) 
    d=CGOLD*e 
  end if 
  u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 ) 
  !! Arrive here with d computed either from parabolic ﬁt, or else from golden
  !! section. 
  fu=func(u) 
  !! This is the one function evaluation per iteration. 
  if (fu <= fx) then 
    !! Now we have to decide what to do with our function evaluation. 
    !! House keeping follows: 
    if (u>= x) then 
      a=x 
    else 
      b=x 
    end if 
    call shft(v,w,x,u) 
    call shft(fv,fw,fx,fu) 
  else 
    if (u < x) then 
      a=u 
    else 
      b=u 
    end if 
    if (fu <= fw .or. w == x) then 
      v=w 
      fv=fw 
      w=u 
      fw=fu 
    else if (fu <= fv .or. v == x .or. v == w) then 
      v=u 
      fv=fu 
    end if 
  end if 
end do !! Done with housekeeping. Back for another iteration. 
call nrerror('brent: exceed maximum iterations') 
CONTAINS 

SUBROUTINE shft(a,b,c,d) 
REAL(SP), INTENT(OUT) :: a 
REAL(SP), INTENT(INOUT) :: b,c 
REAL(SP), INTENT(IN) :: d 
a=b 
b=c 
c=d 
END SUBROUTINE shft 
END FUNCTION brent 

