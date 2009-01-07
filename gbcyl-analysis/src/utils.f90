module utils
  use nrtype, only: dp

  contains 



  !! Author Juho Lintuvuori.
  !! Modified by Jouni Karjalainen to use nrtype dp and intrinsic dot_product
  !!
  !! The angle is defined counterclockwise when looking at the direction of 
  !! (nx, ny, nz).
  !!
  subroutine rotate_vector(x, y, z, nx, ny, nz, angle, xp, yp, zp)
  !
  ! Rotates the vector (x,y,z) into (xp,yp,zp) around axis
  ! (nx,ny,nz) [unit vector of the direction] through angle @p angle
  ! Goldstein: Classical Mechanics 2nd ed., p. 165
  !
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



  subroutine crossp(ax,ay,az,bx,by,bz,cx,cy,cz)
  !
  ! calculates (ax,ay,az) x (bx,by,bz) = cz,cy,cz)
  !
    implicit none
    double precision ax,ay,az,bx,by,bz,cx,cy,cz
  !
    cx = ay * bz - az * by
    cy = az * bx - ax * bz
    cz = ax * by - ay * bx
    return
  end subroutine crossp

end module utils
