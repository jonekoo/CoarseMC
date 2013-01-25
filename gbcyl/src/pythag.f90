FUNCTION pythag_sp(a,b) 
USE nrtype 
IMPLICIT NONE 
REAL(SP),INTENT(IN)::a,b 
REAL(SP)::pythag_sp 
!! Computes (a2 + b2)**(1/2) without destructive underﬂow or overﬂow. 
REAL(SP) :: absa, absb 
absa=abs(a) 
absb=abs(b) 
if(absa>absb)then 
  pythag_sp=absa*sqrt(1.0_sp+(absb/absa)**2) 
else 
  if(absb==0.0)then 
    pythag_sp=0.0 
  else 
    pythag_sp=absb*sqrt(1.0_sp+(absa/absb)**2) 
  endif 
endif 
END FUNCTION pythag_sp 

FUNCTION pythag_dp(a,b) 
USE nrtype 
IMPLICIT NONE 
REAL(DP),INTENT(IN)::a,b 
REAL(DP)::pythag_dp 
REAL(DP)::absa,absb 
absa=abs(a) 
absb=abs(b) 
if(absa>absb)then 
  pythag_dp=absa*sqrt(1.0_dp+(absb/absa)**2) 
else 
  if(absb==0.0)then  
    pythag_dp=0.0 
  else 
    pythag_dp=absb*sqrt(1.0_dp+(absa/absb)**2) 
  endif 
endif 
END FUNCTION pythag_dp 
