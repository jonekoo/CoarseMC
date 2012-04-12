FUNCTION rf_s(x,y,z) 
USE nrtype; USE nrutil, ONLY : assert 
IMPLICIT NONE 
REAL(SP), INTENT(IN) :: x,y,z 
REAL(SP) ::rf_s 
REAL(SP), PARAMETER :: ERRTOL=0.0025_sp !!,TINY=1.5e-38_sp,
REAL(SP), PARAMETER :: BIG=3.0e37_sp,& 
  THIRD=1.0_sp/3.0_sp,&
  C1=1.0_sp/24.0_sp,C2=0.1_sp,C3=3.0_sp/44.0_sp,C4=1.0_sp/14.0_sp 
  !! Computes Carlson?s elliptic integral of the first kind, RF(x,y,z). 
  !! x, y, and z must be nonnegative, and at most one can be zero. TINY must be   
  !! at least 5 times the machine underflow limit, BIG at most one-fifth the
  !! machine overflow limit. 
REAL(SP) ::alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt 
call assert(min(x,y,z) >=0.0, min(x+y,x+z,y+z) >= 5._sp*TINY(1._sp), & 
max(x,y,z) <= BIG, 'rf_s args') 
xt=x
yt=y
zt=z
do
  sqrtx=sqrt(xt) 
  sqrty=sqrt(yt) 
  sqrtz=sqrt(zt)
  alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz 
  xt=0.25_sp*(xt+alamb)
  yt=0.25_sp*(yt+alamb) 
  zt=0.25_sp*(zt+alamb) 
  ave=THIRD*(xt+yt+zt)
  delx=(ave-xt)/ave
  dely=(ave-yt)/ave
  delz=(ave-zt)/ave
  if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit 
end do 
e2=delx*dely-delz**2 
e3=delx*dely*delz 
rf_s=(1.0_sp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
END FUNCTION rf_s  
