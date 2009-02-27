module utils
use nrtype



  contains 



  SUBROUTINE XVEC2(X, Y, Z, NX, NY, NZ, PHI, XP, YP, ZP)
  !
  ! rotates the vector (X,Y,Z) into (XP,YP,ZP) around axis
  ! (NX,NY,NZ) [unit vector of the direction] through angle PHI
  ! Goldstein: Classical Mechanics 2nd ed., p. 165
  !
    IMPLICIT NONE
    REAL(DP), intent(in) :: X, Y, Z
    real(dp), intent(in) :: NX, NY, NZ 
    real(dp), intent(in) :: PHI
    real(dp), intent(out) :: XP, YP, ZP
    REAL(DP) :: DOTP
    DOTP = dot_product((/NX, NY, NZ/), (/X, Y, Z/))
    CALL CROSSP(X, Y, Z, NX, NY, NZ, XP, YP, ZP)
    XP = X * COS(PHI) + NX * DOTP * (1.0 - COS(PHI)) + XP * SIN(PHI)
    YP = Y * COS(PHI) + NY * DOTP * (1.0 - COS(PHI)) + YP * SIN(PHI)
    ZP = Z * COS(PHI) + NZ * DOTP * (1.0 - COS(PHI)) + ZP * SIN(PHI)
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
