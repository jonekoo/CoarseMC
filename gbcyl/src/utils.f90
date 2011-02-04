module utils
use nrtype
use mtmod, only: grnd
USE nr, only: jacobi_sp => jacobi
implicit none

contains 

!> Returns the components of the orientation vector in cylindrical coordinates
!!
  !! @param orientation the orientation vector.
!! @param position the position of the orientation vector.
!! @return the orientation vector in cylindrical coordinates.
!! 
pure function unitvec(orientation, position)
  intrinsic atan2
  real(dp), dimension(3), intent(in) :: orientation
  real(dp), dimension(3), intent(in), optional :: position
  real(dp), dimension(3) :: unitvec
  real(dp) :: uro, utheta, uz
  real(dp) :: nx, ny, nz, theta
  theta = atan2(position(2), position(1))
  nx = 0._dp
  ny = 0._dp
  nz = 1._dp
  call rotate_vector(orientation(1), orientation(2), orientation(3), nx, ny, nz, &
  theta, uro, utheta, uz)
  unitvec = (/uro, utheta, uz/)
end function
  

!> Rotates the vector (x,y,z) into (xp,yp,zp) around axis
!! (nx,ny,nz) [unit vector of the direction] through angle @p angle
!! Goldstein: Classical Mechanics 2nd ed., p. 165
!!
!! Author Juho Lintuvuori.
!! Modified by Jouni Karjalainen to use nrtype dp and intrinsic dot_product
!!
!! The angle is defined counterclockwise when looking at the direction of 
!! (nx, ny, nz).
!!
pure subroutine rotate_vector(x, y, z, nx, ny, nz, angle, xp, yp, zp)
  implicit none
  intrinsic dot_product
  intrinsic cos
  intrinsic sin
  real(dp), intent(in) :: x, y, z, nx, ny, nz, angle
  real(dp), intent(out) :: xp, yp, zp
  real(dp) :: dotpr, cp, sp
  real(dp), dimension(3) :: vec1 
  real(dp), dimension(3) :: vec2 
  vec1 = (/nx, ny, nz/)
  vec2 = (/x, y, z/)
  dotpr = dot_product(vec1, vec2)
  cp = cos(angle)
  sp = sin(angle)
  call crossp(x, y, z, nx, ny, nz, xp, yp, zp)
  xp = x * cp + nx * dotpr * (1.0 - cp) + xp * sp
  yp = y * cp + ny * dotpr * (1.0 - cp) + yp * sp
  zp = z * cp + nz * dotpr * (1.0 - cp) + zp * sp
end subroutine rotate_vector

subroutine nvec(nx, ny, nz)
  !!
  !! Generates a random unit vector (nx, ny, nz). 
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

pure function crossproduct(a, b) result(c)
  !
  ! calculates (ax,ay,az) x (bx,by,bz) = (cz,cy,cz)
  !
  real(dp), dimension(3), intent(in) :: a, b
  real(dp), dimension(3) :: c
  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = a(3) * b(1) - a(1) * b(3)
  c(3) = a(1) * b(2) - a(2) * b(1)
end function

!! Returns a formatting character for the default integer type. To be used when
!! consistent formatting for integer output is needed.
!! 
pure function fmt_char_int() result(format_char)
  character(len = 50) :: format_char
  integer :: r
  character(len = 50) :: w_char
  r = range(r)
  write(w_char, *) 1 + r !! sign + range
  write(format_char, *) 'I' // trim(adjustl(w_char))
  format_char = trim(adjustl(format_char))
end function

!! Returns a formatting character for a double precision real number. To be 
!! used when consistent formatting of real numbers is needed. 
!!
pure function fmt_char_dp() result(format_char)
  integer :: e
  integer :: w
  integer :: d
  real(dp) :: u 
  character(len = 50) :: w_char !! Width of field
  character(len = 50) :: d_char !! Width of decimal field
  character(len = 50) :: e_char !! width of exponent
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

end module
