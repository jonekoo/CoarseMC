module utils
use nrtype
use mtmod, only: grnd
implicit none

contains 

subroutine jacobi(a, n, np, d, v, nrot)
  implicit none
  !Numerical recipes 2nd Edition, Volume 1, page 460
  integer :: n,np,nrot ! np>=n
  real(sp) :: a(np,np), d(np), v(np,np)
  integer, parameter :: nmax=500
  integer :: i,ip,iq,j
  real(sp) :: c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)

  !Alustetaan v yksikkömatriisiksi
  do ip=1,n
     do iq=1,n
        v(ip,iq)=0._sp;
     end do
     v(ip,ip)=1._sp;
  end do
  !Alustetaan b ja d a:n diagonaaliksi
  do ip=1,n
     b(ip)=a(ip,ip);
     d(ip)=b(ip);
     z(ip)=0._sp;
  end do
  nrot=0;
  i=0;
  do
     i=i+1;
     if(i>=50)then
        write(*,*)'liian monta iteraatiota jacobissa!'
        return;
     end if
     sm=0._sp;
     do ip=1,n-1
        do iq=ip+1,n
           sm=sm+abs(a(ip,iq)); !ei diagonaalisten elementtien summa
        end do
     end do
     if(sm==0._sp)return;
     if(i < 4)then
        tresh=0.2_sp*sm/real(n*n, sp); !Ensimmäiset kolme kierrosta
     else
        tresh=0._sp;
     end if
     do ip=1,n-1
        do iq=ip+1,n
           g=100._sp*abs(a(ip,iq));
           !! jos ei-diagonaalinen elementti on riittävän pieni 
           !! jätetään kierto väliin neljännen kierroksen jälkeen
           if((i>4) .and. (abs(d(ip))+g==abs(d(ip))) .and. &
               (abs(d(iq))+g==abs(d(iq))))then
              a(ip,iq)=0._sp;
           else if(abs(a(ip,iq))>tresh)then
              h=d(iq)-d(ip);  
              if(abs(h)+g == abs(h))then
                 t=a(ip,iq)/h;
              else
                 theta=0.5_sp*h/a(ip,iq);
                 t=1._sp/(abs(theta)+sqrt(1._sp+theta*theta));
                 if(theta<0._sp)t=-t;
              end if
              c=1._sp/sqrt(1._sp+t*t);
              s=t*c;
              tau=s/(1._sp+c);
              h=t*a(ip,iq);
              z(ip)=z(ip)-h;
              z(iq)=z(iq)+h;
              d(ip)=d(ip)-h;
              d(iq)=d(iq)+h;
              a(ip,iq)=0._sp;
              do j=1,ip-1
                 g=a(j,ip);
                 h=a(j,iq);
                 a(j,ip)=g-s*(h+g*tau);
                 a(j,iq)=h+s*(g-h*tau);
              end do
              do j=ip+1,iq-1
                 g=a(ip,j);
                 h=a(j,iq);
                 a(ip,j)=g-s*(h+g*tau);
                 a(j,iq)=h+s*(g-h*tau);
              end do
              do j=iq+1,n
                 g=a(ip,j);
                 h=a(iq,j);
                 a(ip,j)=g-s*(h+g*tau);
                 a(iq,j)=h+s*(g-h*tau);
              end do
              do j=1,n
                 g=v(j,ip);
                 h=v(j,iq);
                 v(j,ip)=g-s*(h+g*tau);
                 v(j,iq)=h+s*(g-h*tau);
              end do
              nrot=nrot+1;
           end if
        end do
     end do
     do ip=1,n
        b(ip)=b(ip)+z(ip);
        d(ip)=b(ip);
        z(ip)=0._sp;
     end do
  end do
end subroutine jacobi

!! Returns the components of the orientation vector in cylindrical coordinates
!!
!! :NOTE: It is not possible to call this routine with 
!! call unitvec(particle, particle%ux, particle%uy, particle%uz)
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
  
!! Author Juho Lintuvuori.
!! Modified by Jouni Karjalainen to use nrtype dp and intrinsic dot_product
!!
!! The angle is defined counterclockwise when looking at the direction of 
!! (nx, ny, nz).
!!
pure subroutine rotate_vector(x, y, z, nx, ny, nz, angle, xp, yp, zp)
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

end module utils
