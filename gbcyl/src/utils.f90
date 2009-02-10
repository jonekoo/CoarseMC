module utils
use nrtype

  contains 

  SUBROUTINE XVEC2(X,Y,Z,NX,NY,NZ,PHI,XP,YP,ZP)
  !
  ! rotates the vector (X,Y,Z) into (XP,YP,ZP) around axis
  ! (NX,NY,NZ) [unit vector of the direction] through angle PHI
  ! Goldstein: Classical Mechanics 2nd ed., p. 165
  !
  IMPLICIT NONE
  !
  REAL(DP) X,Y,Z,NX,NY,NZ,PHI,XP,YP,ZP,DP,CP,SP
  !
  ! called functions
  !
  !REAL(DP) DOTP
  !
  DP = dot_product((/NX, NY, NZ/), (/X, Y, Z/))
  CP = COS(PHI)
  SP = SIN(PHI)
   
  CALL CROSSP(X,Y,Z,NX,NY,NZ,XP,YP,ZP)
  XP = X * CP + NX * DP * (1.0 - CP) + XP * SP
  YP = Y * CP + NY * DP * (1.0 - CP) + YP * SP
  ZP = Z * CP + NZ * DP * (1.0 - CP) + ZP * SP
!
!  RETURN
END SUBROUTINE XVEC2



  SUBROUTINE CROSSP(AX,AY,AZ,BX,BY,BZ,CX,CY,CZ)
  !
  ! calculates (AX,AY,AZ) x (BX,BY,BZ) = (CZ,CY,CZ)
  !
    IMPLICIT NONE
    REAL(DP) AX,AY,AZ,BX,BY,BZ,CX,CY,CZ
  !
    CX = AY * BZ - AZ * BY
    CY = AZ * BX - AX * BZ
    CZ = AX * BY - AY * BX
    RETURN
  END SUBROUTINE CROSSP

end module utils
