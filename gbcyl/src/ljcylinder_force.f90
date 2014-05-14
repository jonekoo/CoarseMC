pure function ljcylinder_force(eps, density, sigma, alpha, r, rwall) result(f)
   use nrtype, only: dp
   use ljcylinder, only: EllipticK, EllipticE
   implicit none
   real(dp), intent(in) :: eps, density, sigma, alpha, r, rwall
   real(dp) :: f
   real(dp) :: k
   real(dp), parameter :: Pi = 4._dp * atan(1._dp)
   k = r / rwall

f = (density*Pi*sigma**9*((35 + 5143*k**2 + 27590*k**4 + 27590*k**6 + 5143*k**8 +&
35*k**10 - (160*alpha*(-1 + k**2)**6*(1 + 14*k**2 +&
k**4)*rwall**6)/sigma**6)*EllipticE(k**2) - (-1 + k**2)*(-35 - 3428*k**2 - 15234*k**4&
- 12356*k**6 - 1715*k**8 + (160*alpha*(-1 + k**2)**6*(1 +&
7*k**2)*rwall**6)/sigma**6)*EllipticK(k**2)))/(80.*k*(-1 + k**2)**10*rwall**9)
   f = eps * sigma**3 * f / rwall
end function


pure function EllipticK(k_sqrd) result(Kk)
  use nrtype, only: sp, dp
  use nr, only: rf
  implicit none
  real(dp), intent(in) :: k_sqrd
  real(dp) :: Kk
  real(dp) :: cc, q, rfk, rdk
  cc = 0._dp
  q = 1._dp - k_sqrd
  rfk = real(rf(real(cc, sp), real(q, sp), 1._sp), dp)
  Kk = rfk 
end function 


pure function EllipticE(k_sqrd) result(Ek)
  use nrtype, only: sp, dp
  use nr, only: rd, rf
  real(dp), intent(in) :: k_sqrd
  real(dp) :: Ek
  real(dp) :: cc, q, rfk, rdk
  cc = 0._dp
  q = 1._dp - k_sqrd
  rfk = real(rf(real(cc, sp), real(q, sp), 1._sp), dp)
  rdk = real(rd(real(cc, sp), real(q, sp), 1._sp), dp)
  Ek = rfk - k_sqrd * rdk / 3._dp 
end function 


