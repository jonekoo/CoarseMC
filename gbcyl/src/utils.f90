module utils
use nrtype
use mtmod, only: grnd
implicit none

contains 

subroutine nvec(nx, ny, nz)
  !!
  !! Forms a random unit vector. 
  !!
  !! @see Understanding Mol. Sim. 2nd Ed.  Frenkel, Smit p. 578
  !!
  intrinsic sqrt
  double precision, intent(out) :: nx, ny, nz
  double precision :: l, u1, u2, s
  l = 0.0_dp
  do
     u1 = 1._dp - 2._dp * grnd()
     u2 = 1._dp - 2._dp * grnd()
     l = u1 * u1 + u2 * u2
     if(l <= 1._dp) exit
  end do
  s = 2.0_dp * sqrt(1._dp - l)
  nx = u1 * s
  ny = u2 * s
  nz = 1._dp - 2._dp * l
end subroutine

pure SUBROUTINE XVEC2(X, Y, Z, NX, NY, NZ, PHI, XP, YP, ZP)
  intrinsic cos
  intrinsic sin
  intrinsic dot_product
  !!
  !! rotates the vector (X,Y,Z) into (XP,YP,ZP) around axis
  !! (NX,NY,NZ) [unit vector of the direction] through angle PHI
  !! Goldstein: Classical Mechanics 2nd ed., p. 165
  !!
  REAL(DP), intent(in) :: X, Y, Z
  real(dp), intent(in) :: NX, NY, NZ 
  real(dp), intent(in) :: PHI
  real(dp), intent(out) :: XP, YP, ZP
  REAL(DP) :: DOTP
  DOTP = dot_product((/NX, NY, NZ/), (/X, Y, Z/))
  CALL CROSSP(X, Y, Z, NX, NY, NZ, XP, YP, ZP)
  XP = X * COS(PHI) + NX * DOTP * (1._dp - COS(PHI)) + XP * SIN(PHI)
  YP = Y * COS(PHI) + NY * DOTP * (1._dp - COS(PHI)) + YP * SIN(PHI)
  ZP = Z * COS(PHI) + NZ * DOTP * (1._dp - COS(PHI)) + ZP * SIN(PHI)
END SUBROUTINE

pure SUBROUTINE CROSSP(AX,AY,AZ,BX,BY,BZ,CX,CY,CZ)
  !!
  !! calculates (AX,AY,AZ) x (BX,BY,BZ) = (CZ,CY,CZ)
  !!
  REAL(DP), intent(in) :: AX, AY, AZ, BX, BY, BZ
  real(dp), intent(out) :: CX, CY, CZ
  CX = AY * BZ - AZ * BY
  CY = AZ * BX - AX * BZ
  CZ = AX * BY - AY * BX
END SUBROUTINE

pure function fmt_char_int() result(format_char)
  character(len = 50) :: format_char
  integer :: r
  character(len = 50) :: w_char
  r = range(r)
  write(w_char, *) 1 + r !! sign + range
  write(format_char, *) 'I' // trim(adjustl(w_char))
  format_char = trim(adjustl(format_char))
end function

pure function fmt_char_dp() result(format_char)
  integer :: e
  integer :: w
  integer :: d
  real(dp) :: u 
  character(len = 50) :: w_char
  character(len = 50) :: d_char
  character(len = 50) :: e_char
  character(len = 50) :: format_char
  e = int(log10(real(range(u)))) + 1
  d = precision(u)
  w = 3 + d + 2 + e
  write(w_char, *) w
  write(e_char, *) e
  write(d_char, *) d
  write(format_char, *) 'G' // trim(adjustl(w_char)) // '.' // trim(adjustl(d_char)) // 'E' // trim(adjustl(e_char)) 
  format_char = trim(adjustl(format_char))
end function

end module utils
