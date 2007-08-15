function attwall(xa,rc,sigwall)
use nrtype
use nr, only : rf,rd
implicit none
real(dp), intent(in) :: xa,rc,sigwall
real(dp) :: attwall
real(dp) :: coeff,z,Ek,Kk
real(sp) :: rfk,rdk
real(sp) :: cc,q
z=xa**2
cc=0.0
q=1.0-z
rfk=rf(cc,q,1.0_sp)
rdk=rd(cc,q,1.0_sp)
Ek=rfk-z*rdk/3.0_dp
Kk=rfk
coeff=4.0_dp*atan(1.0)*sigwall/(12*(rc*(1.0_dp-z))**3)
attwall=coeff*((7.0_dp+z)*Ek+4.0_dp*(z-1.0_dp)*Kk)
end function attwall
