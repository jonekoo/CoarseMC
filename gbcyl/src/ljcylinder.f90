module ljcylinder
use nrtype
use nr, only : rf, rd
implicit none

contains

!! Returns the potential energy for a Lennard-Jones (LJ) site with respect to 
!! a cylindrical wall consisting of smoothly and evenly distributed LJ 
!! particles.
!! 
!! @p Kw parameter for adjusting the potential well depth. 
!! @p sigma the parameter setting the contact distance for the first layer in
!! the wall and  the LJ site.
!! @p alpha controls the amount of attraction with respect to repulsion. 
!! alpha=1 is the regular LJ interaction while alpha=0 makes the interaction
!! purely repulsive.
!! @p r the distance of LJ site center from the cylinder axis.
!! @p rwall the radius of the cylinder.
!! 
function ljcylinderpotential(Kw, sigma, alpha, r, rwall) result(potential)
  real(dp), intent(in) :: Kw
  real(dp), intent(in) :: sigma
  real(dp), intent(in) :: alpha
  real(dp), intent(in) :: r
  real(dp), intent(in) :: rwall
  real(dp) :: potential
  potential = Kw * (repwall(r/rwall, rwall/sigma) &
  - alpha * attwall(r/rwall, rwall/sigma))
end function

!! Returns the absolute value of the repulsive part of the interaction
!! between an infinitely thick  wall of a cylindrical cavity consisting of 
!! smoothly and evenly distributed Lennard-Jones particles and a LJ particle.
!!
!! @p xa is the ratio r/rc where r is particle's distance from the center of 
!! the cylinder.
!! @p rc is the radius of the cylinder in reduced LJ units.
!!
function repwall(xa, rc)
  real(dp), intent(in) :: xa, rc !, sigwall
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

!! Returns the absolute value of the attractive part of the interaction
!! between an infinitely thick  wall of a cylindrical cavity consisting of 
!! smoothly and evenly distributed Lennard-Jones particles and a LJ particle.
!!
!! @p xa is the ratio r/rc where r is particle's distance from the center of 
!! the cylinder.
!! @p rc is the radius of the cylinder in reduced LJ units.
!!
function attwall(xa, rc)
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

end module
