!> Returns the absolute value of the repulsive part of the interaction
!! between an infinitely thick  wall of a cylindrical cavity consisting of 
!! smoothly and evenly distributed Lennard-Jones particles and a LJ particle.
!!
!! @param k is the ratio r/rc where r is particle's distance from the center of 
!! the cylinder.
!! @param rc is the radius of the cylinder in reduced LJ units.
!!
pure function repwall2(k, rc)
  use nrtype, only: dp, sp
  use nr, only: rf, rd
  use utils, only: horner
  implicit none
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

  k_sqrd = k**2
  cc = 0._dp
  q = 1._dp - k_sqrd
  rfk = real(rf(real(cc, sp), real(q, sp), 1._sp), dp)
  rdk = real(rd(real(cc, sp), real(q, sp), 1._sp), dp)
  Ek = rfk - k_sqrd * rdk / 3._dp !! Complete elliptic integrals of the second
  Kk = rfk                        !! and first kind

  C_E9 = horner(as, k_sqrd)
  C_K9 = horner(bs, k_sqrd)

  I_10 = 2._dp / (9._dp * ((rc) * (1._dp - k_sqrd))**9) * &
         (C_E9 * Ek + C_K9 * Kk) 

  repwall2 = 63._dp * pi / 64._dp * I_10
end function



!> Returns the absolute value of the attractive part of the interaction
!! between an infinitely thick  wall of a cylindrical cavity consisting of 
!! smoothly and evenly distributed Lennard-Jones particles and a LJ particle.
!!
!! @param k is the ratio r/rc where r is particle's distance from the center of 
!! the cylinder.
!! @param rc is the radius of the cylinder in reduced LJ units.
!!
pure function attwall2(k, rc)
  use nrtype, only: dp, sp
  use nr, only: rf, rd
  implicit none
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

