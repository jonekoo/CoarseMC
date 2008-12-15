!! Calculates the absolute value of the repulsive part of the interaction
!! between a smooth infinitely thick cylinder shape wall and a Lennard-Jones
!! particle inside it. 
!! x=r/rc where r=particles distance from the center of the cylinder,
!! rc=radius of the cylinder
!! sigwall=contact distance between the particle and the wall

function repwall(xa,rc,sigwall)
use nrtype
use nr, only : rf,rd
implicit none
real(dp), intent(in) :: xa,rc,sigwall
real(dp) :: repwall
real(dp), dimension(5) :: as,bs
real(dp) :: pa,pb,z
real(sp) :: rfk,rdk
real(dp) :: Ek,Kk,coeff
real(sp) :: cc,q
  interface 
  function horner(a,n,x) 
    use nrtype
    real(dp), dimension(:), intent(in) :: a
    integer, intent(in) :: n 
    real(dp), intent(in) :: x
    real(dp) :: horner
  end function horner
  end interface

z=xa**2
cc=0.0
q=1.0-z
rfk=rf(cc,q,1.0_sp)
rdk=rd(cc,q,1.0_sp)
Ek=rfk-z*rdk/3.0_dp   !! Complete elliptic integrals of the second and
Kk=rfk                      !! first kind
as(1:5)=(/35.0_dp,4052.0_dp,16434.0_dp,11156.0_dp,1091.0_dp/)
bs(1:5)=(/1400.0_dp,6304.0_dp,-1200.0_dp,-5728.0_dp,-776.0_dp/)
pa=horner(as,5,z)
pb=horner(bs,5,z)
coeff=4.0_dp*atan(1.0)*(sigwall**12)/(128.0_dp*45.0_dp*((rc)*(1.0_dp-z))**9)
repwall=coeff*(pa*Ek+pb*Kk)
end function repwall
