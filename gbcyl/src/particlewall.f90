!! Jatkuvasti ja tasaisesti jakautuneista LJ-partikkeleista
!! koostuvan seinän ja LJ-partikkelin väliseen vuorovaikutukseen
!! liittyviä parametreja ja aliohjelmia.  

module particlewall
  use particle, only: particledat 
  use nrtype
  use class_poly_box
  use class_parameterizer
  use class_parameter_writer
  implicit none
  private

  public :: initptwall
  public :: particlewall_potential
  public :: repwall
  public :: attwall
  public :: gbwall
  public :: particlewall_write_parameters
    
  !! @see Micheletti et. al. Journal of Chemical Physics 123, 224705 for 
  !! definitions of these parameters.
  !!
  real(dp), save :: alpha_A_ = 1._dp
  real(dp), save :: alpha_B_ = 1._dp
  real(dp), save :: K_w_ = 8._dp
  
  !! The distance of LJ interaction sites from the Gay-Berne particle center
  !! along the unique axis.
  real(dp), save :: LJdist = 1.7_dp 
  logical, save :: is_uniform_alignment = .false.

  interface initptwall
    module procedure initptwall_parameterizer
  end interface

  interface particlewall_potential
    module procedure gbwall
  end interface

contains 

  subroutine initptwall_parameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call get_parameter(reader, 'Kw', K_w_)
    call get_parameter(reader, 'alpha_A', alpha_A_) 
    call get_parameter(reader, 'alpha_B', alpha_B_) 
    call get_parameter(reader, 'Kw', K_w_)
    call get_parameter(reader, 'LJ_dist', LJdist) 
    call get_parameter(reader, 'is_uniform_alignment', is_uniform_alignment)
  end subroutine

  subroutine particlewall_write_parameters(writer)
    type(parameter_writer), intent(in) :: writer
    call write_comment(writer, 'Particle-wall interaction parameters')
    call write_parameter(writer, 'alpha_A', alpha_A_) 
    call write_parameter(writer, 'alpha_B', alpha_B_) 
    call write_parameter(writer, 'Kw', K_w_)
    call write_parameter(writer, 'LJ_dist', LJdist) 
    call write_parameter(writer, 'is_uniform_alignment', is_uniform_alignment)
  end subroutine

  subroutine gbwall(gbparticle, simbox, Eptwall,ovrlp)
    implicit none
    type(particledat), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: Eptwall
    logical, intent(out) :: ovrlp
    real(dp), parameter :: sig=1._dp
    real(dp) :: rsiteA, rsiteB
    real(dp) :: repulA, attracA, repulB, attracB, fu, Rc
    real(dp) :: rAperR, rBperR
    ovrlp = .false.
    Eptwall = 0._dp
    Rc = get_x(simbox)/2._dp !! :TODO: Remove this dependency by making a cylinder wall object.
    call rArB(gbparticle, rsiteA, rsiteB)
    if(rsiteA >= Rc .or. rsiteB >= Rc) then
      ovrlp = .true.
      return;
    end if
    rAperR  = rsiteA / Rc
    rBperR  = rsiteB / Rc
    attracA = attwall(rAperR, Rc, sig)
    repulA  = repwall(rAperR, Rc, sig)
    attracB = attwall(rBperR, Rc, sig)
    repulB  = repwall(rBperR, Rc, sig)
    if (is_uniform_alignment) then
      fu = angular(gbparticle)  
    else 
      fu = 1._dp
    end if
    Eptwall = fu * K_w_ * ((repulA - alpha_A_ * attracA) + (repulB - alpha_B_ &
      * attracB))     
  end subroutine gbwall

  !! Returns the distances of the interaction sites from the axis of the 
  !! cylinder.
  !! 
  !! @p particle the particle of which interaction sites are calcula 
  !! @p r_a the distance of site A from the cylinder center
  !! @p rb the distance of site B from the cylinder center
  subroutine rarb(particle, r_a, r_b)
    implicit none
    intrinsic sqrt
    type(particledat), intent(in) :: particle
    real(dp), intent(out) :: r_a, r_b
    real(dp) :: x_a, y_a, x_b, y_b
    x_a=particle%x+LJdist*particle%ux
    y_a=particle%y+LJdist*particle%uy
    x_b=particle%x-LJdist*particle%ux
    y_b=particle%y-LJdist*particle%uy
    r_a=sqrt(x_a**2 + y_a**2)
    r_b=sqrt(x_b**2 + y_b**2)
  end subroutine



  real(dp) function angular(particle)
    implicit none
    type(particledat), intent(in) :: particle
    angular=(particle%uz)**2
  end function angular



!! Calculates the absolute value of the repulsive part of the interaction
!! between a smooth infinitely thick cylinder shape wall and a Lennard-Jones
!! particle inside it. 
!! x=r/rc where r=particles distance from the center of the cylinder,
!! rc=radius of the cylinder
!! sigwall=contact distance between the particle and the wall
!!
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
  coeff = 4._dp * atan(1._dp) * (sigwall**12) / (128._dp * 45._dp * ((rc) * &
  (1._dp - z))**9)
  repwall = coeff * (pa * Ek + pb * Kk)
end function repwall



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



end module particlewall
