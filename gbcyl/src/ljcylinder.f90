module ljcylinder
use nrtype
use nr, only : rf, rd
implicit none


interface 
  pure function ljcylinder_force(eps, density, sigma, alpha, r, rwall) result(f)
    use nrtype, only: dp
    real(dp), intent(in) :: eps, density, sigma, alpha, r, rwall
    real(dp) :: f
  end function
end interface


interface
pure function EllipticE(k_sqrd) result(Ek)
  use nrtype, only: dp
  use nr, only: rd, rf
  real(dp), intent(in) :: k_sqrd
  real(dp) :: Ek
end function
end interface

interface
pure function EllipticK(k_sqrd) result(Kk)
  use nrtype, only: dp
  use nr, only: rf
  implicit none
  real(dp), intent(in) :: k_sqrd
  real(dp) :: Kk
end function
end interface

contains

!! Returns the potential energy for a Lennard-Jones (LJ) site with respect to 
!! a cylindrical wall consisting of smoothly and evenly distributed LJ 
!! particles.
!! 
!! @p eps the parameter for LJ-particle - wall well-depth
!! @p density the wall density
!! @p sigma the parameter setting the contact distance for the first layer in
!! the wall and  the LJ site.
!! @p alpha controls the amount of attraction with respect to repulsion. 
!! alpha=1 is the regular LJ interaction while alpha=0 makes the interaction
!! purely repulsive.
!! @p r the distance of LJ site center from the cylinder axis.
!! @p rwall the radius of the cylinder.
!! 
pure function ljcylinderpotential(eps, density, sigma, alpha, r, rwall) result(potential)
  real(dp), intent(in) :: eps 
  real(dp), intent(in) :: density 
  real(dp), intent(in) :: sigma
  real(dp), intent(in) :: alpha
  real(dp), intent(in) :: r
  real(dp), intent(in) :: rwall
  real(dp) :: potential
  potential = eps * density * sigma**3 * (repwall2(r/rwall, rwall/sigma) &
  - alpha * attwall2(r/rwall, rwall/sigma))
end function


!! Returns the absolute value of the repulsive part of the interaction
!! between an infinitely thick  wall of a cylindrical cavity consisting of 
!! smoothly and evenly distributed Lennard-Jones particles and a LJ particle.
!!
!! @p xa is the ratio r/rc where r is particle's distance from the center of 
!! the cylinder.
!! @p rc is the radius of the cylinder in reduced LJ units.
!!
pure function repwall(xa, rc)
  real(dp), intent(in) :: xa, rc !, sigwall
  real(dp) :: repwall
  real(dp), dimension(5) :: as, bs
  real(dp) :: pa, pb, z
  real(dp) :: rfk, rdk
  real(dp) :: Ek, Kk, coeff
  real(dp) :: cc, q
  interface 
  pure function horner(a, n, x) 
    use nrtype
    real(dp), dimension(:), intent(in) :: a
    integer, intent(in) :: n 
    real(dp), intent(in) :: x
    real(dp) :: horner
  end function horner
  end interface
  z = xa**2
  cc = 0._dp
  q = 1._dp - z
  rfk = real(rf(real(cc, sp), real(q, sp), 1._sp), dp)
  rdk = real(rd(real(cc, sp), real(q, sp), 1._sp), dp)
  Ek = rfk - z * rdk / 3._dp   !! Complete elliptic integrals of the second and
  Kk = rfk                     !! first kind
  !! as are the coefficients for the polynomial in equation (11) and
  !! bs are the coefficients for the polynomial in equation (12).
  !! pa and pb are the corresponding polynomial values. 
  !! They have been multiplied by 315 and divided by two just to make things 
  !! look a bit cleaner.  
  !! Coefficient of the highest order term is first on the array as required by
  !! the horner function. 
  as(1:5) = (/35._dp, 4052._dp, 16434._dp, 11156._dp, 1091._dp/)
  bs(1:5) = (/1400._dp, 6304._dp, -1200._dp, -5728._dp, -776._dp/)
  pa = horner(as, 5, z)
  pb = horner(bs, 5, z)
  coeff = 4._dp * atan(1._dp) / (128._dp * 45._dp * ((rc) * &
  (1._dp - z))**9)
  repwall = coeff * (pa * Ek + pb * Kk)
end function repwall

!! Returns the absolute value of the repulsive part of the interaction
!! between an infinitely thick  wall of a cylindrical cavity consisting of 
!! smoothly and evenly distributed Lennard-Jones particles and a LJ particle.
!!
!! @p k is the ratio r/rc where r is particle's distance from the center of 
!! the cylinder.
!! @p rc is the radius of the cylinder in reduced LJ units.
!!
pure function repwall2(k, rc)
  real(dp), intent(in) :: k, rc !, sigwall
  real(dp) :: repwall2
  real(dp) :: rfk, rdk
  real(dp) :: Ek, Kk
  real(dp) :: cc, q, k_sqrd
  real(dp), parameter :: pi = 4._dp * atan(1._dp)

  !! as are the coefficients for the polynomial in equation (11) and
  !! bs are the coefficients for the polynomial in equation (12).
  !! Coefficient of the highest order term is first on the array as required by
  !! the horner function. 
  real(dp), parameter :: as(5) = (/2._dp / 9._dp, 8104._dp / 315._dp, &
  3652._dp / 35._dp, 22312._dp / 315._dp, 2182._dp / 315._dp/)
  real(dp), parameter :: bs(5) = (/80._dp / 9._dp, 12608._dp / 315._dp, &
  -160._dp / 21._dp, -11456._dp / 315._dp, -1552._dp / 315._dp /)
  !! C_E9 and C_K9 are the corresponding polynomial values. 
  real(dp) :: C_E9, C_K9
  real(dp) :: I_10 !! As defined by equation (16) with m = 10.

  interface 
  pure function horner(a, n, x) 
    use nrtype
    real(dp), dimension(:), intent(in) :: a
    integer, intent(in) :: n 
    real(dp), intent(in) :: x
    real(dp) :: horner
  end function horner
  end interface

  k_sqrd = k**2
  cc = 0._dp
  q = 1._dp - k_sqrd
  rfk = real(rf(real(cc, sp), real(q, sp), 1._sp), dp)
  rdk = real(rd(real(cc, sp), real(q, sp), 1._sp), dp)
  Ek = rfk - k_sqrd * rdk / 3._dp !! Complete elliptic integrals of the second
  Kk = rfk                        !! and first kind

  C_E9 = horner(as, 5, k_sqrd)
  C_K9 = horner(bs, 5, k_sqrd)

  I_10 = 2._dp / (9._dp * ((rc) * (1._dp - k_sqrd))**9) * &
         (C_E9 * Ek + C_K9 * Kk) 

  repwall2 = 63._dp * pi / 64._dp * I_10
end function


!! Returns the absolute value of the attractive part of the interaction
!! between an infinitely thick  wall of a cylindrical cavity consisting of 
!! smoothly and evenly distributed Lennard-Jones particles and a LJ particle.
!!
!! @p xa is the ratio r/rc where r is particle's distance from the center of 
!! the cylinder.
!! @p rc is the radius of the cylinder in reduced LJ units.
!!
pure function attwall(xa, rc)
  real(dp), intent(in) :: xa, rc !, sigwall
  real(dp) :: attwall
  real(dp) :: coeff, z, Ek, Kk
  real(dp) :: rfk, rdk
  real(dp) :: cc, q
  z = xa**2
  cc = 0._dp
  q = 1._dp - z
  rfk = real(rf(real(cc, sp), real(q, sp), 1._sp), dp)
  rdk = real(rd(real(cc, sp), real(q, sp), 1._sp), dp)
  Ek = rfk - z * rdk / 3._dp
  Kk = rfk
  coeff = 4._dp * atan(1._dp) / (12._dp * (rc * (1._dp - z))**3)
  attwall = coeff * ((7._dp + z) * Ek + 4._dp * (z - 1._dp) * Kk)
end function attwall


!! Returns the absolute value of the attractive part of the interaction
!! between an infinitely thick  wall of a cylindrical cavity consisting of 
!! smoothly and evenly distributed Lennard-Jones particles and a LJ particle.
!!
!! @p k is the ratio r/rc where r is particle's distance from the center of 
!! the cylinder.
!! @p rc is the radius of the cylinder in reduced LJ units.
!!
pure function attwall2(k, rc)
  real(dp), intent(in) :: k, rc !, sigwall
  real(dp) :: attwall2
  real(dp) :: k_sqrd, Ek, Kk
  real(dp) :: rfk, rdk
  real(dp) :: cc, q
  real(dp) :: I_4         !! As defined by eq. (16) with m = 4.
  real(dp) :: C_E3, C_K3  !! Polynomials of eqs. (13) and (14), respectively. 
  real(dp), parameter :: pi = 4._dp * atan(1._dp)
  k_sqrd = k**2
  cc = 0._dp
  q = 1._dp - k_sqrd
  rfk = real(rf(real(cc, sp), real(q, sp), 1._sp), dp)
  rdk = real(rd(real(cc, sp), real(q, sp), 1._sp), dp)
  Ek = rfk - k_sqrd * rdk / 3._dp
  Kk = rfk

  C_E3 = (14._dp + 2._dp * k_sqrd) / 3._dp 
  C_K3 = 8._dp * (k_sqrd - 1._dp) / 3._dp 

  I_4 = 2._dp / (3._dp * (rc * (1._dp - k_sqrd))**3) * (C_E3 * Ek + C_K3 * Kk)

  attwall2 = (3._dp * pi / 2._dp) * I_4
end function

end module
