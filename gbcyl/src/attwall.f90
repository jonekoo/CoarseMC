function attwall(xa, rc, sigwall)
use nrtype
use nr, only : rf, rd
implicit none
real(dp), intent(in) :: xa, rc, sigwall
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
coeff = 4._dp * atan(1._dp) * sigwall / (12._dp * (rc * (1._dp - z))**3)
attwall = coeff * ((7._dp + z) * Ek + 4._dp * (z - 1._dp) * Kk)
end function attwall
