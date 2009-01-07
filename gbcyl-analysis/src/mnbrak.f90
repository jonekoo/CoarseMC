SUBROUTINE mnbrak(ax, bx, cx, fa, fb, fc, func) 
USE nrtype; USE nrutil, ONLY : swap 
IMPLICIT NONE 
REAL(SP), INTENT(INOUT) :: ax,bx 
REAL(SP), INTENT(OUT) :: cx,fa,fb,fc 
INTERFACE 
  FUNCTION func(x) 
  USE nrtype 
  IMPLICIT NONE 
  REAL(SP), INTENT(IN) :: x 
  REAL(SP) :: func 
  END FUNCTION func 
END INTERFACE 
REAL(SP), PARAMETER :: GOLD=1.618034_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp 
  !! Given a function func, and given distinct initial points ax and bx, this 
  !! routine searches in the downhill direction (deﬁned by the function as 
  !! evaluated at the initial points) and returns new points ax, bx, cx that 
  !! bracket a minimum of the function. Also returned are the function values 
  !! at the three points, fa, fb, and fc. 
  !! Parameters: GOLD is the default ratio by which successive intervals are 
  !! magniﬁed; GLIMIT is the maximum magniﬁcation allowed for a parabolic-ﬁt 
  !! step. 
REAL(SP) ::fu,q,r,u,ulim 
fa=func(ax) 
fb=func(bx) 
if (fb > fa) then 
  !! Switch roles of a and b so that we can go downhill in the direction 
  !! from a to b. 
  call swap(ax, bx) 
  call swap(fa, fb) 
end if 
cx=bx+GOLD*(bx-ax) !! First guess for c. 
fc=func(cx) 
do 
  !! Do-while-loop: Keep returning here until we bracket. 
  if (fb < fc) RETURN 
  !! Compute u by parabolic extrapolation from a, b, c. TINY is used to 
  !! prevent any possible division by zero. 
  r=(bx-ax)*(fb-fc) 
  q=(bx-cx)*(fb-fa) 
  u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_sp*sign(max(abs(q-r),TINY),q-r)) 
  ulim=bx+GLIMIT*(cx-bx) 
  !! We won’t go farther than this. Test various possibilities: 
  if ((bx-u)*(u-cx) > 0.0) then 
    !! Parabolic u is between b and c: try it. 
    fu=func(u) 
    if (fu <fc) then 
      !! Got a minimum between b and c. 
      ax=bx 
      fa=fb 
      bx=u 
      fb=fu 
      RETURN 
    else if (fu > fb) then 
      !! Got a minimum between a and u. 
      cx=u 
      fc=fu 
      RETURN 
    end if 
    u=cx+GOLD*(cx-bx) !! Parabolic ﬁt was no use. Use default magniﬁcation. 
    fu=func(u) 
  else if ((cx-u)*(u-ulim) > 0.0) then 
    !! Parabolic ﬁt is between c and its allowed limit. 
    fu=func(u) 
    if (fu <fc) then 
      bx=cx 
      cx=u 
      u=cx+GOLD*(cx-bx) 
      call shft(fb,fc,fu,func(u)) 
    end if 
  else if ((u-ulim)*(ulim-cx) >= 0.0) then 
    !! Limit parabolic u to maximum allowed value. 
    u=ulim 
    fu=func(u) 
  else 
    !! Reject parabolic u, use default magniﬁcation. 
    u=cx+GOLD*(cx-bx) 
    fu=func(u) 
  end if 
  call shft(ax,bx,cx,u) 
  call shft(fa,fb,fc,fu) !! Eliminate oldest point and continue. 
end do 
CONTAINS 

SUBROUTINE shft(a,b,c,d) 
REAL(SP), INTENT(OUT) :: a 
REAL(SP), INTENT(INOUT) :: b,c 
REAL(SP), INTENT(IN) :: d 
a=b 
b=c 
c=d 
END SUBROUTINE shft 
END SUBROUTINE mnbrak 

