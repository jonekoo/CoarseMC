!! Jatkuvasti ja tasaisesti jakautuneista LJ-partikkeleista
!! koostuvan sein�n ja LJ-partikkelin v�liseen vuorovaikutukseen
!! liittyvi� parametreja ja aliohjelmia.  

module particlewall
  use particle, only: particledat 
  use nrtype
  use class_poly_box
  use class_parameterizer
  use class_parameter_writer
  implicit none
  private

  public :: new_sitewallpotential
  public :: new_gbwallpotential
  public :: initptwall
  public :: particlewall_potential
  public :: repwall
  public :: attwall
  public :: particlewall_writeparameters
  public :: energy
  public :: sitewallpotential
  public :: gbwallpotential
    
  !! @see Micheletti et. al. Journal of Chemical Physics 123, 224705 for 
  !! definitions of these parameters.
  !!
  real(dp), save :: alphaA = 1._dp
  real(dp), save :: alphaB = 1._dp
  real(dp), save :: Kw = 8._dp
  real(dp), save :: sig=1._dp
  
  !! The distance of LJ interaction sites from the Gay-Berne particle center
  !! along the unique axis.
  real(dp), save :: LJdist = 1.7_dp 
  logical, save :: isuniformalignment = .false.

!! LJ-GB interaction consists of interactions of the two sites. Basicly
!! That could be in a different module which then would use this module. 
!! Keep in mind that the Xe-wall interaction may be of the same form as 
!! the GB-wall interaction. Things that characterize a typical Lennard-Jones
!! interaction are the contact distance sigma and the potential well-depth 
!! parameter epsilon. What is added here is basicly only the alpha parameter. 
!! In addition to that the strength parameter Kw includes both the well depth
!! and the density of the wall plus some additional constants.  

type sitewallpotential
  private
  real(dp) :: alpha = 1._dp    ! relative strength of attraction 
  real(dp) :: sigma = 1._dp    ! LJ contact distance 
  real(dp) :: Kw = 8._dp       ! interaction strength 
end type

type gbwallpotential
  private
  type(sitewallpotential) :: sitea
  type(sitewallpotential) :: siteb
  real(dp) :: separation = 1.7_dp
  logical :: isuniformalignment = .true.
end type

interface initptwall
  module procedure initptwallparameterizer, initptwallregular
end interface

interface particlewall_potential
  module procedure gbwall
end interface

interface energy
  module procedure sitewallpotential_energy
end interface

contains 

function new_sitewallpotential(Kw, sigma, alpha) result(potential)
  real(dp), intent(in), optional :: Kw
  real(dp), intent(in), optional :: sigma
  real(dp), intent(in), optional :: alpha
  type(sitewallpotential) :: potential
  if(present(Kw)) potential%Kw = Kw
  if(present(sigma)) potential%sigma = sigma
  if(present(alpha)) potential%alpha = alpha
end function

function new_gbwallpotential(sitea, siteb, separation, isuniformalignment) &
result(potential)
  type(sitewallpotential), intent(in), optional :: sitea
  type(sitewallpotential), intent(in), optional :: siteb
  real(dp), intent(in), optional :: separation
  logical, intent(in), optional :: isuniformalignment
  type(gbwallpotential) :: potential
  if(present(sitea)) potential%sitea = sitea
  if(present(siteb)) potential%siteb = siteb
  if(present(separation)) potential%separation = separation
  if(present(isuniformalignment)) potential%isuniformalignment = &
  isuniformalignment 
end function

subroutine initptwallregular(Kwin, alphaAin, alphaBin, LJdistin, isunifin, sigwallin)
  real(dp), intent(in), optional :: Kwin
  real(dp), intent(in), optional :: alphaAin
  real(dp), intent(in), optional :: alphaBin
  real(dp), intent(in), optional :: LJdistin
  logical, intent(in), optional :: isunifin
  real(dp), intent(in), optional :: sigwallin
  if(present(Kwin)) Kw = Kwin
  if(present(alphaAin)) alphaA = alphaAin
  if(present(alphaBin)) alphaB = alphaBin
  if(present(LJdistin)) LJdist = LJdistin
  if(present(isunifin)) isuniformalignment = isunifin
  if(present(sigwallin)) sig = sigwallin
end subroutine

  subroutine initptwallparameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'Kw', Kw)
    call getparameter(reader, 'alpha_A', alphaA) 
    call getparameter(reader, 'alpha_B', alphaB) 
    call getparameter(reader, 'LJ_dist', LJdist) 
    call getparameter(reader, 'is_uniform_alignment', isuniformalignment)
    call getparameter(reader, 'sigwall', sig)
  end subroutine

  subroutine particlewall_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'Particle-wall interaction parameters')
    call writeparameter(writer, 'alpha_A', alphaA) 
    call writeparameter(writer, 'alpha_B', alphaB) 
    call writeparameter(writer, 'Kw', Kw)
    call writeparameter(writer, 'LJ_dist', LJdist) 
    call writeparameter(writer, 'is_uniform_alignment', isuniformalignment)
  end subroutine

  subroutine gbwall(gbparticle, simbox, Eptwall,ovrlp)
    implicit none
    type(particledat), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: Eptwall
    logical, intent(out) :: ovrlp
    real(dp) :: rsiteA, rsiteB
    real(dp) :: repulA, attracA, repulB, attracB, fu, Rc
    real(dp) :: rAperR, rBperR
    ovrlp = .false.
    Eptwall = 0._dp
    !! :TODO: Remove this dependency by making a cylinder wall object.
    Rc = getx(simbox)/2._dp 
    call rArB(gbparticle, rsiteA, rsiteB)
    if(rsiteA >= Rc .or. rsiteB >= Rc) then
      ovrlp = .true.
      return
    end if
    rAperR  = rsiteA / Rc
    rBperR  = rsiteB / Rc
    attracA = attwall(rAperR, Rc, sig)
    repulA  = repwall(rAperR, Rc, sig)
    attracB = attwall(rBperR, Rc, sig)
    repulB  = repwall(rBperR, Rc, sig)
    if (isuniformalignment) then
      fu = angular(gbparticle)  
    else 
      fu = 1._dp
    end if
    Eptwall = fu * Kw * ((repulA - alphaA * attracA) + (repulB - alphaB &
      * attracB))     
  end subroutine gbwall



subroutine sitewallpotential_energy(sitepotential, r, rwall, energy, overlap)
  type(sitewallpotential), intent(in) :: sitepotential
  real(dp), intent(in) :: r
  real(dp), intent(in) :: rwall
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  energy = 0._dp
  overlap = .false.
  if (r >= rwall) then
    overlap = .true.
    return
  end if
  energy = sitepotential%Kw * (repwall(r/rwall, rwall, sitepotential%sigma) &
  - sitepotential%alpha * attwall(r/rwall, rwall, sitepotential%sigma))
end subroutine

  !! Returns the distances of the interaction sites from the axis of the 
  !! cylinder.
  !! 
  !! @p particle the particle of which interaction sites are calcula 
  !! @p ra the distance of site A from the cylinder center
  !! @p rb the distance of site B from the cylinder center
  subroutine rarb(particle, ra, rb)
    implicit none
    intrinsic sqrt
    type(particledat), intent(in) :: particle
    real(dp), intent(out) :: ra, rb
    real(dp) :: xa, ya, xb, yb
    xa=particle%x+LJdist*particle%ux
    ya=particle%y+LJdist*particle%uy
    xb=particle%x-LJdist*particle%ux
    yb=particle%y-LJdist*particle%uy
    ra=sqrt(xa**2 + ya**2)
    rb=sqrt(xb**2 + yb**2)
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
