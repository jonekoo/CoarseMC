!> Module responsible for defining the computation of the (uniaxial)
!! Gay-Berne potential for the interaction of two ellipsoidal molecules.
!!
!! @see J. G. Gay and B. J. Berne. J. Chem. Phys., 74(6):3316, 1981.
!! @see M. A. Bates and G. R. Luckhurst. J. Chem. Phys.,
!!      110(14):7087, 1999.
!!
module gayberne
use nrtype, only: sp, dp
use class_parameterizer
use class_parameter_writer
implicit none
private

public :: gayberne_init
public :: gb_potential
public :: sigma
public :: getsigma0
public :: getkappasigma
public :: gb_writeparameters
public :: d_potential
public :: gb_R
public :: g_potential
public :: gb_force
public :: gb_epsilon
 
!> Ratio of contact distances for end-to-end and side-by-side
!! configurations of two particles.
!! @see Luckhurst & et.al J.Chem.Phys, Vol. 110, No. 14
real(dp), save :: kappasigma   = 4.4_dp   !! = sige/sigs

!> Ratio of well-depths for side-by-side and end-to-end
!! configurations of two particles.
real(dp), save :: kappaepsilon = 20._dp !! = epss/epse  

!> This parameter is for fine-tuning how much the side-by-side configuration
!! of two particles is favored as compared to the end-to-end configuration.  
real(dp), save :: mu            = 1._dp

!> This parameter is for fine-tuning how much parallel alignment of two
!! particles is favored as compared to other configurations
!! (e.g. cross-configuration).
real(dp), save :: nu            = 1._dp

!> The contact distance (where potential is zero) at cross-configuration
!! of two particles.
real(dp), save :: sigma0       = 1._dp

!> The well-depth at the cross-configuration of two particles. 
real(dp), save :: epsilon0     = 1._dp


!> Parameters below are derived from parameters above. These are saved
!! to increase performance.
real(dp), save :: chiepsilon
real(dp), save :: chisigma
real(dp), save :: chisigmasquared

!> Defines the GB distance when the particles overlap in gb_potential.
real(dp), save :: gb_hardcore = 0.6_dp

!> Defines the procedure to be used for computing the well-depth
!! function epsilon. This allows common special cases of the
!! parameterization to be treated more efficiently than a single
!! general function.
procedure(gb_eps), pointer :: gb_epsilon

!> Initializes the module.
interface gayberne_init
  module procedure initparameterizer, initold
end interface

!> The interface for the well-depth function in the GB potential.
interface
   pure function gb_eps(ui, uj, urij)
     use nrtype
     implicit none
     real(dp), intent(in) :: ui(3), uj(3), urij(3)
     real(dp) :: gb_eps
   end function gb_eps
end interface

contains

!> Initializes the module using a parameterizer object.
!! 
!! @param[in] reader the parameterizer object which is responsible for getting the
!! parameters for this module from some source, e.g. file. 
!!
subroutine initparameterizer(reader)
  type(parameterizer), intent(in) :: reader
  call getparameter(reader, 'gb_kappa_sigma', kappasigma)
  call getparameter(reader, 'gb_kappa_epsilon', kappaepsilon)
  call getparameter(reader, 'gb_mu', mu)
  call getparameter(reader, 'gb_nu', nu)
  call getparameter(reader, 'gb_sigma_0', sigma0)
  call getparameter(reader, 'gb_epsilon_0', epsilon0)
  call getparameter(reader, 'gb_hardcore', gb_hardcore)
  call init_common()
end subroutine initparameterizer

subroutine init_common()
  if (abs(mu - 1) < tiny(mu) .and. abs(nu - 1) < tiny(nu)) then
     write(*, *) '# gb_nu = 1 and gb_nu = 1. Using gb_epsilon_mu1nu1'
     gb_epsilon => gb_epsilon_mu1nu1
  else
     gb_epsilon => gb_epsilon_full
  end if
  chiepsilon = &
       (kappaepsilon**(1._dp / mu) - 1._dp) / &
       (kappaepsilon**(1._dp / mu)+ 1._dp)
  chisigma = (kappasigma * kappasigma - 1._dp) / &
       (kappasigma * kappasigma + 1._dp)
  chisigmasquared = chisigma**2
end subroutine init_common


!> Initializes the module for potential calculation. 
!! 
!! @see M. A. Bates and G. R. Luckhurst, JCP 110(14), 7078, 1999 for a
!! detailed discussion.
!!
!! @param kappasigmain sets the axis ratio sigma_ee/sigma_ss of the
!!        ellipsoidal molecule.
!! @param kappaepsilonin sets the ratio of well depths in side-by-side
!!        and end-to-end configurations epsilon_ss/epsilon_ee.
!! @param muin adjusts how much side by side configuration is favored.
!! @param nuin parameter for adjusting how much parallel alignment in
!!        favored.
!! @param sigma0in sets the contact distance for two ellipsoids in a
!!        cross configuration. For a one-component Gay-Berne liquid this
!!        can be set to 1.
!! @param epsilon0in sets the well depth for the potential. For a one
!!        component Gay-Berne liquid this can be set to 1.
!! @param hardcore can be given to set a hard core to the potential. In effect
!!        this defines the gb_R at which two particles overlap.
!! 
subroutine initold(kappasigmain, kappaepsilonin, muin, nuin, sigma0in, &
     epsilon0in, hardcore)
  real(dp), intent(in) :: kappasigmain
  real(dp), intent(in) :: kappaepsilonin
  real(dp), intent(in) :: muin
  real(dp), intent(in) :: nuin
  real(dp), intent(in) :: sigma0in
  real(dp), intent(in) :: epsilon0in
  real(dp), intent(in), optional :: hardcore
  kappasigma = kappasigmain
  kappaepsilon = kappaepsilonin
  mu = muin
  nu = nuin
  sigma0 = sigma0in
  epsilon0 = epsilon0in
  if (present(hardcore)) gb_hardcore = hardcore
  call init_common()
end subroutine initold


!> Write the parameters of this module to the output unit and format
!! defined by @p writer.
subroutine gb_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call writecomment(writer, 'Gay-Berne potential parameters')
  call writeparameter(writer, 'gb_kappa_sigma', kappasigma)
  call writeparameter(writer, 'gb_kappa_epsilon', kappaepsilon)
  call writeparameter(writer, 'gb_mu', mu)
  call writeparameter(writer, 'gb_nu', nu)
  call writeparameter(writer, 'gb_sigma_0', sigma0)
  call writeparameter(writer, 'gb_epsilon_0', epsilon0)    
  call writeparameter(writer, 'gb_hardcore', gb_hardcore)
end subroutine gb_writeparameters


!> Calculates the Gay-Berne potential for two particles. If @p overlap
!! is true, the potential is zero.
!!
!! @param ui,uj unit orientation vectors of the particles.
!! @param rij vector from the center of particle i to particle j.
!! @param energy the interaction energy.
!! @param overlap true if the potential has a hard core and the particles
!!        are too close to each other. 
!!
pure subroutine gb_potential(ui, uj, rij, energy, overlap)
  implicit none
  real(dp), dimension(3), intent(in) :: rij, ui, uj
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  real(dp) :: rgb, gb6
  real(dp), dimension(3) :: urij
  energy = 0.0
  overlap = .false.
  rgb = gb_R(ui, uj, rij)
  if(rgb < gb_hardcore) then
     overlap = .true.
  else
     gb6 = rgb**(-6)
     energy = gb6 * (gb6 - 1._dp)
     urij = rij / sqrt(dot_product(rij, rij))
     energy = 4._dp * gb_epsilon(ui, uj, urij) * energy
  end if
end subroutine gb_potential


!> Returns the Gay-Berne potential for two particles. Overlap of
!! particles is not considered. 
!!
!! @param ui,uj the unit orientation vectors for the particles.
!! @param rij the vector from the center of particle i to the center of
!!        particle j.
!! 
real(dp) pure function potentialf(ui, uj, rij) result(energy)
  real(dp), dimension(3), intent(in) :: rij, ui, uj
  real(dp) :: rgb, gb6
  real(dp), dimension(3) :: urij
  rgb = gb_R(ui, uj, rij)
  gb6 = rgb**(-6)
  energy = gb6 * (gb6 - 1._dp)
  urij = rij / sqrt(dot_product(rij, rij))
  energy = 4._dp * gb_epsilon(ui, uj, urij) * energy
end function potentialf


!> Calculates the derivative of the potential with respect to the
!! @p alpha coordinate of the vector @p rij.
real(dp) pure function d_potential(ui, uj, rij, alpha)
  real(dp), dimension(3), intent(in) :: rij, ui, uj
  integer, intent(in) :: alpha
  real(dp) :: rijabs
  real(dp), dimension(3) :: urij
  real(dp) :: energy
  logical :: overlap
  rijabs = sqrt(dot_product(rij, rij))
  urij = rij/rijabs
  call gb_potential(ui, uj, rij, energy, overlap)
  d_potential = energy / gb_epsilon(ui, uj, urij) * mu * & 
       rd_anisotropic2(ui, uj, urij, chiepsilon, alpha) / rijabs + &
       4._dp * gb_epsilon(ui, uj, urij) * (6._dp * gb_R(ui, uj, rij)**(-7) - &
       12._dp * gb_R(ui, uj, rij)**(-13)) / &
       sigma0 * (rij(alpha) / rijabs + sigma0 / 2._dp / &
       sqrt(anisotropic(ui, uj, urij, chisigma))**3 * &
       rd_anisotropic2(ui, uj, urij, chisigma, alpha) / rijabs)
end function

real(dp) pure function rd_anisotropic2(ui, uj, urij, chi, alpha)
  real(dp), dimension(3), intent(in) :: ui, uj, urij
  real(dp), intent(in) :: chi
  integer, intent(in) :: alpha
  rd_anisotropic2 = 2 * urij(alpha) * (1 - anisotropic(ui, uj, urij, chi)) - &
       chi / (1 - chi**2 * dot_product(ui,uj)**2) * &
       (tijalpha(ui, ui, urij, alpha) + tijalpha(uj, uj, urij, alpha) - &
       2 * chi * dot_product(ui,uj) * tijalpha(ui, uj, urij, alpha))
end function

real(dp) pure function tijalpha(ui, uj, urij, alpha)
  real(dp), dimension(3), intent(in) :: ui, uj, urij
  integer, intent(in) :: alpha
  tijalpha = dot_product(uj, urij) * ui(alpha) + dot_product(ui, urij) * uj(alpha)
end function

!> The anisotropic distance function of the GB potential.
!!
!! @param ui,uj the unit orientation vectors of the two particles.
!! @param rij the distance vector from the center of particle i to the center
!!        of particle j.
!!
real(dp) pure function gb_R(ui, uj, rij)
  implicit none
  real(dp), dimension(3), intent(in) :: ui
  real(dp), dimension(3), intent(in) :: uj
  real(dp), dimension(3), intent(in) :: rij
  real(dp) :: r
  r = sqrt(dot_product(rij, rij))
  gb_R = (r - sigma(ui, uj, rij / r) + sigma0) / sigma0
end function gb_R

real(dp) pure function anisotropic(ui, uj, urij, chi)
  real(dp), dimension(3), intent(in) :: ui
  real(dp), dimension(3), intent(in) :: uj
  real(dp), dimension(3), intent(in) :: urij
  real(dp), intent(in) :: chi
  real(dp) :: idj
  real(dp) :: ids
  real(dp) :: jds
  idj = dot_product(ui, uj) 
  ids = dot_product(ui, urij)
  jds = dot_product(uj, urij)
  anisotropic = 1 - chi * (ids**2 + jds**2 - 2 * chi * ids * jds * idj) / &
       (1._dp - chi**2 * idj**2) 
end function anisotropic

!> Calculates the derivative of function anisotropic with respect to the 
!! @p alpha component of the vector rij = @p urij * rijabs
real(dp) pure function d_anisotropic(ui, uj, urij, rijabs, chi, alpha)
  real(dp), dimension(3), intent(in) :: ui
  real(dp), dimension(3), intent(in) :: uj
  real(dp), dimension(3), intent(in) :: urij
  real(dp), intent(in) :: rijabs 
  real(dp), intent(in) :: chi
  integer, intent(in) :: alpha
  real(dp) :: idj
  real(dp) :: ids
  real(dp) :: jds
  real(dp), dimension(3) :: idt, jdt
  real(dp), dimension(3, 3) :: t  
  integer :: i, j
  idj = dot_product(ui, uj) 
  ids = dot_product(ui, urij)
  jds = dot_product(uj, urij)
  forall(i=1:3, j=1:3) t(i, j) = -urij(i)*urij(j)
  forall(i=1:3) t(i,i) = t(i,i) + 1._dp
  t(:,:) = t(:,:) / rijabs
  idt = matmul(t, ui)
  jdt = matmul(t, uj)
  d_anisotropic = -2 * chi * (ids * idt(alpha) + jds * idt(alpha) - &
       chi * idj * idt(alpha) * jds - chi * idj * ids * jdt(alpha)) / & 
       (1 - chi**2 * idj**2)
end function d_anisotropic


!> Calculates the derivative of function anisotropic with respect to the 
!! @p alpha component of the vector rij = @p urij * rijabs
real(dp) pure function rd_anisotropic(ui, uj, urij, chi, alpha)
  real(dp), dimension(3), intent(in) :: ui
  real(dp), dimension(3), intent(in) :: uj
  real(dp), dimension(3), intent(in) :: urij
  real(dp), intent(in) :: chi
  integer, intent(in) :: alpha
  real(dp) :: idj
  real(dp) :: ids
  real(dp) :: jds
  idj = dot_product(ui, uj) 
  ids = dot_product(ui, urij)
  jds = dot_product(uj, urij)
  rd_anisotropic = -2 * chi / (1 - chi**2 * idj**2) * (ids * ui(alpha) - &
       ids**2 * urij(alpha) + jds * uj(alpha) - jds**2 * urij(alpha) - &
       chi * idj * (ids * ui(alpha) + jds * uj(alpha) - 2 * ids * jds * &
       urij(alpha)))
end function rd_anisotropic


!> The anisotropic shape function of the GB potential.
!! 
!! @param ui,uj the unit orientations of the two particles.
!! @param urij the unit vector pointing from the center of particle i to
!!        the center of particle j.
!! 
real(dp) pure function sigma(ui, uj, urij)
  real(dp), dimension(3), intent(in) :: ui
  real(dp), dimension(3), intent(in) :: uj
  real(dp), dimension(3), intent(in) :: urij
  sigma = sigma0 / sqrt(anisotropic(ui, uj, urij, chisigma))
end function sigma

real(dp) pure function sigmahelp(ids, jds, idj)
  real(dp), intent(in) :: ids, jds, idj
  real(dp) :: idssq, jdssq, idjsq
  idssq = ids * ids
  jdssq = jds * jds
  idjsq = idj * idj
  sigmahelp = 1 - chisigma * (idssq + jdssq - 2 * chisigma * ids * jds * idj)&
       / (1 - chisigmasquared * idjsq)
  sigmahelp = sigma0 / sqrt(sigmahelp)
end function sigmahelp

real(dp) pure function gb_epsilon_full(ui, uj, urij)
  real(dp), dimension(3), intent(in) :: ui
  real(dp), dimension(3), intent(in) :: uj
  real(dp), dimension(3), intent(in) :: urij
  real(dp) :: idj
  idj = dot_product(ui, uj) 
  gb_epsilon_full = epsilon0 * ep(idj)**nu * anisotropic(ui, uj, urij, chiepsilon)**mu
end function gb_epsilon_full


real(dp) pure function gb_epsilon_mu1nu1(ui, uj, urij)
  real(dp), dimension(3), intent(in) :: ui
  real(dp), dimension(3), intent(in) :: uj
  real(dp), dimension(3), intent(in) :: urij
  real(dp) :: idj
  idj = dot_product(ui, uj) 
  gb_epsilon_mu1nu1 = epsilon0 * ep(idj) * anisotropic(ui, uj, urij, chiepsilon)
end function gb_epsilon_mu1nu1


real(dp) pure function ep(idj)
  real(dp), intent(in) :: idj
  ep = 1._dp / sqrt(1._dp - chisigmasquared * idj**2)
end function ep

real(dp) pure function epp(ids, jds, idj)
  real(dp),intent(in) :: ids, jds, idj
  epp = 1 - chiepsilon * (ids**2 + jds**2 - 2 * chiepsilon * ids * jds * idj)&
       / (1 - chiepsilon**2 * idj**2)
end function epp

!> Returns the contact distance for two GB particles in a cross-configuration.
real(dp) function getsigma0()
  getsigma0 = sigma0
end function getsigma0

!> Returns the ratio of contact distances in the end-to-end and
!! side-by-side configurations. 
real(dp) function getkappasigma()
  getkappasigma = kappasigma
end function getkappasigma

!! below here the "new" functions for force calculation

real(dp) pure function d_anisotropic_ids(ui, uj, urij, chi)
  real(dp), dimension(3), intent(in) :: ui, uj, urij
  real(dp), intent(in) :: chi
  d_anisotropic_ids = -2._dp*chi/(1._dp-chi**2*dot_product(ui, uj)**2)*&
       (dot_product(ui, urij)-chi*dot_product(uj, urij)*dot_product(ui, uj))
end function d_anisotropic_ids

real(dp) pure function d_sigma_anisotropic(aniso)
  real(dp), intent(in) :: aniso
  d_sigma_anisotropic = -0.5_dp*sigma0*sqrt(aniso)**(-3)
end function d_sigma_anisotropic

real(dp) pure function d_gb_R_sigma()
  d_gb_R_sigma = -1._dp/sigma0
end function d_gb_R_sigma

real(dp) pure function d_potential_gb_R(gb_eps, gb_r)
  real(dp), intent(in) :: gb_eps
  real(dp), intent(in) :: gb_r
  d_potential_gb_R = 4._dp*gb_eps*(6._dp*gb_r**(-7)-12._dp*gb_r**(-13))
end function d_potential_gb_R

real(dp) pure function d_potential_epp(epp, pot)
  real(dp), intent(in) :: epp
  real(dp), intent(in) :: pot
  d_potential_epp =  pot*mu/epp
end function d_potential_epp

pure function g_ids(ui, urij, rijabs)
  real(dp), dimension(3), intent(in) :: ui, urij
  real(dp), intent(in) :: rijabs
  real(dp), dimension(3) :: g_ids 
  g_ids = (ui - dot_product(urij, ui)*urij)/rijabs
end function g_ids

!> Returns the gradient of the Gay-Berne potential (with respect to
!! @p rij).
!!
!! @param ui,uj the unit orientation vectors of the two particles.
!! @param rij the vector from the center of particle i to the center of
!!        particle j.
pure function g_potential(ui, uj, rij)
  real(dp), dimension(3), intent(in) :: ui, uj, rij
  real(dp) :: rijabs
  real(dp), dimension(3) :: g_potential
  real(dp) :: d_gb_R_rijabs 
  real(dp), dimension(3) :: urij
  real(dp) :: ids, jds, idj
  real(dp), dimension(3) :: g_rijabs
  real(dp), parameter :: d_epp_anisotropic = 1._dp
  d_gb_R_rijabs = 1.0 / sigma0
  rijabs = sqrt(dot_product(rij, rij))
  urij = rij / rijabs
  g_rijabs = urij
  ids = dot_product(ui, urij)
  jds = dot_product(uj, urij)
  idj = dot_product(ui, uj)
  g_potential = d_potential_gb_R(gb_epsilon(ui, uj, urij), gb_R(ui, uj, rij)) &
       * d_gb_R_rijabs * g_rijabs + &
       d_potential_epp(epp(ids, jds, idj), potentialf(ui, uj, rij)) * &
       d_epp_anisotropic * (d_anisotropic_ids(ui, uj, urij, chiepsilon) * &
       g_ids(ui, urij, rijabs) + d_anisotropic_ids(uj, ui, urij, chiepsilon) &
       * g_ids(uj, urij, rijabs)) + &
       d_potential_gb_R(gb_epsilon(ui, uj, urij), gb_R(ui, uj, rij)) * &
       d_gb_R_sigma() * &
       d_sigma_anisotropic(anisotropic(ui, uj, urij, chisigma))*( &
       d_anisotropic_ids(ui, uj, urij, chisigma)*g_ids(ui, urij, rijabs) + &
       d_anisotropic_ids(uj, ui, urij, chisigma)*g_ids(uj, urij, rijabs))
  end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure function grad_anisotropic(ui, uj, rij, chi)
  real(dp) :: grad_anisotropic(3)
  real(dp), intent(in) :: ui(3), uj(3), rij(3), chi
  real(dp) :: mat_u(3, 3)
  integer :: k, l
  forall(k=1:3, l=1:3) mat_u(k, l) = ui(k) * ui(l) + uj(k) * uj(l) - &
    2 * chi * dot_product(ui, uj) * ui(k) * uj(l)
  grad_anisotropic = matmul((mat_u + transpose(mat_u)), rij) / &
       dot_product(rij, rij) + dot_product(matmul(rij, mat_u), rij) * (-2) / &
       dot_product(rij, rij)**2 * rij
  grad_anisotropic = grad_anisotropic * (-chi) / &
       (1 - chi**2 * dot_product(ui, uj)**2)
end function


!> Returns the Gay-Berne force between two particles with unit
!! orientation vectors @p ui,uj and center-to-center vector @p rij.
pure function gb_force(ui, uj, rij)
  real(dp), intent(in) :: ui(3), uj(3), rij(3)
  real(dp) :: gb_force(3)
  real(dp) :: urij(3)
  real(dp) :: e
  real(dp) :: grad_e(3)
  real(dp) :: grad_gb_sigma(3)
  real(dp) :: grad_gb_R(3)
  real(dp) :: r

  urij = rij / sqrt(dot_product(rij, rij))

  !! Potential well depth and its gradient:
  e = gb_epsilon(ui, uj, urij)
  grad_e = epsilon0 * ep(dot_product(ui, uj)) * &
       grad_anisotropic(ui, uj, rij, chiepsilon)

  !! The gradients of the Gay-Berne contact distance and distance
  grad_gb_sigma = -0.5 * sigma0 * &
       sqrt(anisotropic(ui, uj, urij, chisigma))**(-3) * &
       grad_anisotropic(ui, uj, rij, chisigma)
  grad_gb_R = (urij - grad_gb_sigma) / sigma0
  r = gb_R(ui, uj, rij)

  gb_force = -4 * (grad_e * (r**(-12) - r**(-6)) + &
    e * (6 * r**(-7) - 12 * r**(-13)) * grad_gb_R)

end function

end module


