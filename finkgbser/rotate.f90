module rotate
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
  DOUBLE PRECISION X,Y,Z,NX,NY,NZ,PHI,XP,YP,ZP,DP,CP,SP
  !
  ! called functions
  !
  !DOUBLE PRECISION DOTP
  !
  DP = DOTP(NX,NY,NZ,X,Y,Z)
  CP = COS(PHI)
  SP = SIN(PHI)
   
  CALL CROSSP(X,Y,Z,NX,NY,NZ,XP,YP,ZP)
  XP = X * CP + NX * DP * (1.0 - CP) + XP * SP
  YP = Y * CP + NY * DP * (1.0 - CP) + YP * SP
  ZP = Z * CP + NZ * DP * (1.0 - CP) + ZP * SP
!
!  RETURN
END SUBROUTINE XVEC2


  FUNCTION DOTP(AX,AY,AZ,BX,BY,BZ) result(dotprod)
  !
  ! calculates (AX,AY,AZ) . (BX,BY,BZ) = DOTP
  !
    IMPLICIT NONE
    real(dp),intent(in) ::  AX,AY,AZ,BX,BY,BZ
    real(dp) :: dotprod 
  !
    dotprod = AX * BX + AY * BY + AZ * BZ
    
    RETURN
  END FUNCTION DOTP



  SUBROUTINE CROSSP(AX,AY,AZ,BX,BY,BZ,CX,CY,CZ)
  !
  ! calculates (AX,AY,AZ) x (BX,BY,BZ) = CZ,CY,CZ)
  !
    IMPLICIT NONE
    DOUBLE PRECISION AX,AY,AZ,BX,BY,BZ,CX,CY,CZ
  !
    CX = AY * BZ - AZ * BY
    CY = AZ * BX - AX * BZ
    CZ = AX * BY - AY * BX
    RETURN
  END SUBROUTINE CROSSP



  subroutine nvec(genrand,nx,ny,nz)
    ! Aliohjelma satunnaisen yksikkövektorin muodostamista varten
    ! Understanding Mol. Sim. 2nd Ed.  Frenkel, Smit
    ! s.  578
    implicit none
    double precision,intent(out) :: nx,ny,nz
    interface 
      real function genrand() 
      end function genrand
    end interface
    double precision :: l,u1,u2,s
    l=0.0;
    do
       u1=1.0-2.0*genrand();
       u2=1.0-2.0*genrand();
       l=u1*u1+u2*u2;
       if(l<=1.0)exit;
    end do
    s=2.0*sqrt(1.0-l);
    nx=u1*s;
    ny=u2*s;
    nz=1.0-2.0*l;
  end subroutine nvec


  subroutine rotangle(genrand,dthetamax,theta)
    implicit none
    double precision,intent(in) :: dthetamax
    double precision,intent(out) :: theta
    interface
      real function genrand() 
      end function genrand
    end interface
    theta=(2.0*genrand()-1.0)*dthetamax;
  end subroutine rotangle
  
end module rotate
