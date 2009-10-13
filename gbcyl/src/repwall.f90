!! Calculates the absolute value of the repulsive part of the interaction
!! between a smooth infinitely thick cylinder shape wall and a Lennard-Jones
!! particle inside it. 
!! x=r/rc where r=particles distance from the center of the cylinder,
!! rc=radius of the cylinder
!! sigwall=contact distance between the particle and the wall

function repwall(xa, rc, sigwall)
use nrtype
use nr, only : rf, rd
implicit none
real(dp), intent(in) :: xa, rc, sigwall
real(dp) :: repwall
real(dp), dimension(5) :: as, bs
real(dp) :: pa, pb, z
real(dp) :: rfk, rdk
real(dp) :: Ek, Kk, coeff
real(dp) :: cc, q
interface 
function horner(a, n, x) 
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
  as(1:5) = (/35._dp, 4052._dp, 16434._dp, 11156._dp, 1091._dp/)
  bs(1:5) = (/1400._dp, 6304._dp, -1200._dp, -5728._dp, -776._dp/)
  pa = horner(as, 5, z)
  pb = horner(bs, 5, z)
  coeff = 4._dp * atan(1._dp) * (sigwall**12) / (128._dp * 45._dp * ((rc) * &
  (1._dp - z))**9)
  repwall = coeff * (pa * Ek + pb * Kk)
end function repwall
