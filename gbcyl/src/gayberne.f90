module gayberne
  use nrtype, only: sp, dp
  use class_parameterizer
  use class_parameter_writer
  implicit none
  private

  !! :TODO: Remove dependencies on parameter writing and reading to the pair 
  !! :TODO: potential module. One should be able to use this module without 
  !! :TODO: being dependent on any parameter reading solution. 

  public :: gayberne_init
  public :: potential
  public :: sigma
  public :: getsigma0
  public :: getkappasigma
  public :: gb_epsilon
  public :: gb_writeparameters
  public :: d_potential
  public :: gb_R
  public :: g_potential
  public :: chisigma, chiepsilon
  public :: gblj_potential
 
  !! Parameters of the Gay-Berne potential. Documentation:
  !! @see Luckhurst & et.al J.Chem.Phys, Vol. 110, No. 14
  !!
  real(dp), save :: kappasigma   = 4.4_dp   !! = sige/sigs
  real(dp), save :: kappaepsilon = 20._dp !! = epss/epse  
  real(dp), save :: mu            = 1._dp
  real(dp), save :: nu            = 1._dp
  real(dp), save :: sigma0       = 1._dp
  real(dp), save :: epsilon0     = 1._dp


  !! Parameters below are derived from parameters above.
  real(dp), save :: chiepsilon
  real(dp), save :: chisigma
  real(dp), save :: chisigmasquared

  !! Defines when the particles "overlap"
  real(dp), parameter :: hardcore = 0.6_dp

  !! Parameters for the GB-LJ interaction:
  real(dp), save :: gblj_epsilon_0 = 1._dp
  real(dp), save :: gblj_sigma_0 = 1._dp

  interface potential
    module procedure potentials
  end interface

  interface gayberne_init
    module procedure initparameterizer, initold
  end interface

  contains

  !! Initializes the module using a parameterizer object.
  !! 
  !! @p reader the parameterizer object which is responsible for getting the
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
    call getparameter(reader, 'gblj_epsilon_0', gblj_epsilon_0)
    call getparameter(reader, 'gblj_sigma_0', gblj_sigma_0)
    chiepsilon = &
      (kappaepsilon**(1._dp / mu) - 1._dp) / &
      (kappaepsilon**(1._dp / mu)+ 1._dp)
    chisigma = (kappasigma * kappasigma - 1._dp) / &
      (kappasigma * kappasigma + 1._dp)
    chisigmasquared = chisigma**2
  end subroutine

  !! Initializes the module for potential calculation. 
  !! 
  !! @see M. A. Bates and G. R. Luckhurst, JCP 110(14), 7078, 1999 for a 
  !! detailed discussion.
  !!
  !! @p kappasigmain sets the axis ratio sigma_ee/sigma_ss of the 
  !! ellipsoidal molecule.
  !! @p kappaepsilonin sets the ratio of well depths in side-by-side and 
  !! end-to-end configurations epsilon_ss/epsilon_ee.
  !! @p muin adjusts how much side by side configuration is favored.
  !! @p nuin parameter for adjusting how much parallel alignment in favored. 
  !! @p sigma0in sets the contact distance for two ellipsoids in a cross 
  !! configuration. For a one-component Gay-Berne liquid this can be set to 1.
  !! @p epsilon0in sets the well depth for the potential. For a one component
  !! Gay-Berne liquid this can be set to 1.
  !! 
  subroutine initold(kappasigmain, kappaepsilonin, muin, nuin, sigma0in, epsilon0in)
    implicit none
    real(dp), intent(in) :: kappasigmain
    real(dp), intent(in) :: kappaepsilonin
    real(dp), intent(in) :: muin
    real(dp), intent(in) :: nuin
    real(dp), intent(in) :: sigma0in
    real(dp), intent(in) :: epsilon0in
    kappasigma = kappasigmain
    kappaepsilon = kappaepsilonin
    mu = muin
    nu = nuin
    sigma0 = sigma0in
    epsilon0 = epsilon0in
    chiepsilon = &
      (kappaepsilon**(1._dp / mu) - 1._dp) / &
      (kappaepsilon**(1._dp / mu)+ 1._dp)
    chisigma = (kappasigma * kappasigma - 1._dp) / &
      (kappasigma * kappasigma + 1._dp)
    chisigmasquared = chisigma**2
  end subroutine

  subroutine gb_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'Gay-Berne potential parameters')
    call writeparameter(writer, 'gb_kappa_sigma', kappasigma)
    call writeparameter(writer, 'gb_kappa_epsilon', kappaepsilon)
    call writeparameter(writer, 'gb_mu', mu)
    call writeparameter(writer, 'gb_nu', nu)
    call writeparameter(writer, 'gb_sigma_0', sigma0)
    call writeparameter(writer, 'gb_epsilon_0', epsilon0)    
    call writeparameter(writer, 'gblj_epsilon_0', gblj_epsilon_0)
    call writeparameter(writer, 'gblj_sigma_0', gblj_sigma_0)
  end subroutine

  !! Calculates the Gay-Berne potential for two particles 
  !! particlei and particlej. If rgb=rij-sig+sig0<=0.6 
  !! routine returns with ovrlp=.true. 
  !!
  !! :TODO: check if ovrlp is really needed and in which conditions.
  !! :TODO: make hardcore adjustable
  !!
  pure subroutine potentials(ui, uj, rij, gbV, ovrlp)
    implicit none
    real(dp), dimension(3), intent(in) :: rij, ui, uj
    real(dp), intent(out) :: gbV
    logical, intent(out) :: ovrlp
    real(dp) :: rgb, gb6
    real(dp), dimension(3) :: urij
    gbV = 0._dp
    ovrlp = .false.
    rgb = gb_R(ui, uj, rij)
    if(rgb < hardcore) then
      ovrlp = .true.
    else
      gb6 = rgb**(-6)
      gbV = gb6 * (gb6 - 1._dp)
      urij = rij / sqrt(dot_product(rij, rij))
      gbV = 4._dp * gb_epsilon(ui, uj, urij) * gbV
    end if
  end subroutine

  !! Calculates the Gay-Berne potential for two particles 
  !! particlei and particlej. If rgb=rij-sig+sig0<=0.6 
  !! routine returns with ovrlp=.true. 
  !!
  !! :TODO: check if ovrlp is really needed and in which conditions.
  !! :TODO: make hardcore adjustable
  !!
  real(dp) pure function potentialf(ui, uj, rij) result(gbV)
    implicit none
    real(dp), dimension(3), intent(in) :: rij, ui, uj
    real(dp) :: rgb, gb6
    real(dp), dimension(3) :: urij
    rgb = gb_R(ui, uj, rij)
    gb6 = rgb**(-6)
    gbV = gb6 * (gb6 - 1._dp)
    urij = rij / sqrt(dot_product(rij, rij))
    gbV = 4._dp * gb_epsilon(ui, uj, urij) * gbV
  end function

  !! Calculates the derivative of the potential with respect to the @p alpha 
  !! coordinate of the vector @p rij.
  real(dp) pure function d_potential(ui, uj, rij, alpha)
    real(dp), dimension(3), intent(in) :: rij, ui, uj
    integer, intent(in) :: alpha
    real(dp) :: rijabs
    real(dp), dimension(3) :: urij
    real(dp) :: energy
    logical :: overlap
    rijabs = sqrt(dot_product(rij, rij))
    urij = rij/rijabs
    call potential(ui, uj, rij, energy, overlap)
    d_potential = energy/gb_epsilon(ui, uj, urij)*mu* & !d_anisotropic(ui, uj, urij, rijabs, chiepsilon, alpha)+&
      rd_anisotropic2(ui, uj, urij, chiepsilon, alpha)/rijabs+ &
      4._dp*gb_epsilon(ui, uj, urij)*(6._dp*gb_R(ui, uj, rij)**(-7) - &
      12._dp*gb_R(ui, uj, rij)**(-13))/&
      sigma0*(rij(alpha)/rijabs+sigma0/2._dp/&
      sqrt(anisotropic(ui, uj, urij, chisigma))**3*&
      rd_anisotropic2(ui, uj, urij, chisigma, alpha)/rijabs)
      !d_anisotropic(ui, uj, urij, rijabs, chisigma, alpha))
  end function

real(dp) pure function rd_anisotropic2(ui, uj, urij, chi, alpha)
  real(dp), dimension(3), intent(in) :: ui, uj, urij
  real(dp), intent(in) :: chi
  integer, intent(in) :: alpha
  rd_anisotropic2 = 2._dp*urij(alpha)*(1._dp-anisotropic(ui, uj, urij, chi))-&
    chi/(1._dp-chi**2*dot_product(ui,uj)**2)*&
    (tijalpha(ui, ui, urij, alpha)+tijalpha(uj, uj, urij, alpha)-&
    2._dp*chi*dot_product(ui,uj)*tijalpha(ui, uj, urij, alpha))
end function

real(dp) pure function tijalpha(ui, uj, urij, alpha)
  real(dp), dimension(3), intent(in) :: ui, uj, urij
  integer, intent(in) :: alpha
  tijalpha = dot_product(uj, urij)*ui(alpha) + dot_product(ui, urij)*uj(alpha)
end function

  real(dp) pure function gb_R(ui, uj, rij)
    implicit none
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: rij
    real(dp) :: r
    r = sqrt(dot_product(rij, rij))
    gb_R = (r - sigma(ui, uj, rij/r) + sigma0) / sigma0
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
    anisotropic = 1._dp - chi * (ids**2 + jds**2 - 2._dp * chi * &
      ids * jds * idj) / (1._dp - chi**2 * idj**2) 
  end function

  !! Calculates the derivative of function anisotropic with respect to the 
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
    t(:,:) = t(:,:)/rijabs
    idt = matmul(t, ui)
    jdt = matmul(t, uj)
    d_anisotropic = -2._dp*chi*(ids*idt(alpha)+jds*idt(alpha)-&
      chi*idj*idt(alpha)*jds-chi*idj*ids*jdt(alpha))/(1._dp-chi**2*idj**2)
  end function

  !! Calculates the derivative of function anisotropic with respect to the 
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
    rd_anisotropic = -2._dp*chi/(1._dp-chi**2*idj**2)*(ids*ui(alpha)-&
      ids**2*urij(alpha)+jds*uj(alpha)-jds**2*urij(alpha)-&
      chi*idj*(ids*ui(alpha)+jds*uj(alpha)-2._dp*ids*jds*urij(alpha)))
  end function

  real(dp) pure function sigma(ui, uj, urij)
    implicit none
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: urij
    !real(dp) :: idj
    !real(dp) :: ids
    !real(dp) :: jds
    !idj = dot_product(ui, uj) 
    !ids = dot_product(ui, urij)
    !jds = dot_product(uj, urij)
    !sigma = sigmahelp(ids, jds, idj)
    sigma = sigma0/sqrt(anisotropic(ui, uj, urij, chisigma))
  end function sigma

  !! :TODO: change this to take unit vectors as arguments. 
  real(dp) pure function sigmahelp(ids, jds, idj)
  implicit none
  intrinsic sqrt
  real(dp), intent(in) :: ids, jds, idj
  real(dp) :: idssq, jdssq, idjsq
    idssq = ids*ids
    jdssq = jds*jds
    idjsq = idj*idj
    sigmahelp = 1._dp - chisigma * (idssq + jdssq - 2._dp * chisigma * &
      ids * jds * idj) / (1._dp - chisigmasquared * idjsq)
    sigmahelp = sigma0 / sqrt(sigmahelp)
  end function sigmahelp

  real(dp) pure function gb_epsilon(ui, uj, urij)
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: urij
    real(dp) :: idj
    !real(dp) :: ids
    !real(dp) :: jds
    idj = dot_product(ui, uj) 
    !ids = dot_product(ui, urij)
    !jds = dot_product(uj, urij)
    !gb_epsilon = epsilon0 * ep(idj)**nu * epp(ids, jds, idj)**mu
    gb_epsilon = epsilon0 * ep(idj)**nu * anisotropic(ui, uj, urij, chiepsilon)**mu
  end function


  !! :TODO: change this to take unit vectors as arguments. 
  real(dp) pure function ep(idj)
  implicit none
  intrinsic sqrt
  real(dp), intent(in) :: idj
    ep = 1._dp / sqrt(1._dp - chisigmasquared * idj**2)
  end function ep

  !! :TODO: change this to take unit vectors as arguments.   
  real(dp) pure function epp(ids, jds, idj)
  implicit none
  real(dp),intent(in) :: ids, jds, idj
    epp = 1._dp - chiepsilon * (ids**2 + jds**2 - 2._dp * chiepsilon * &
      ids * jds * idj) / (1._dp - chiepsilon**2 * idj**2)
  end function epp

  function getsigma0()
  real(dp) :: getsigma0
    getsigma0 = sigma0
  end function 

  function getkappasigma()
  real(dp) :: getkappasigma
    getkappasigma = kappasigma
  end function 

  !! below here the new functions for force calculation

  real(dp) pure function d_anisotropic_ids(ui, uj, urij, chi)
    real(dp), dimension(3), intent(in) :: ui, uj, urij
    real(dp), intent(in) :: chi
    d_anisotropic_ids = -2._dp*chi/(1._dp-chi**2*dot_product(ui, uj)**2)*&
      (dot_product(ui, urij)-chi*dot_product(uj, urij)*dot_product(ui, uj))
  end function 

  real(dp) pure function d_sigma_anisotropic(aniso)
    real(dp), intent(in) :: aniso
    d_sigma_anisotropic = -0.5_dp*sigma0*sqrt(aniso)**(-3)
  end function

  real(dp) pure function d_gb_R_sigma()
    d_gb_R_sigma = -1._dp/sigma0
  end function

  real(dp) pure function d_potential_gb_R(gb_eps, gb_r)
    real(dp), intent(in) :: gb_eps
    real(dp), intent(in) :: gb_r
    d_potential_gb_R = 4._dp*gb_eps*(6._dp*gb_r**(-7)-12._dp*gb_r**(-13))
  end function

  real(dp) pure function d_potential_epp(epp, pot)
    real(dp), intent(in) :: epp
    real(dp), intent(in) :: pot
    d_potential_epp =  pot*mu/epp
  end function

  pure function g_ids(ui, urij, rijabs)
    real(dp), dimension(3), intent(in) :: ui, urij
    real(dp), intent(in) :: rijabs
    real(dp), dimension(3) :: g_ids 
    g_ids = (ui - dot_product(urij, ui)*urij)/rijabs
  end function

  pure function g_potential(ui, uj, rij)
    real(dp), dimension(3), intent(in) :: ui, uj, rij
    real(dp) :: rijabs
    real(dp), dimension(3) :: g_potential
    real(dp) :: d_gb_R_rijabs 
    real(dp), dimension(3) :: urij
    real(dp) :: ids, jds, idj
    real(dp), dimension(3) :: g_rijabs
    real(dp), parameter :: d_epp_anisotropic = 1._dp
    d_gb_R_rijabs = 1._dp/sigma0
    rijabs = sqrt(dot_product(rij, rij))
    urij = rij/rijabs
    g_rijabs = urij
    ids = dot_product(ui, urij)
    jds = dot_product(uj, urij)
    idj = dot_product(ui, uj)
    g_potential = d_potential_gb_R(gb_epsilon(ui, uj, urij), gb_R(ui, uj, rij))*&
      d_gb_R_rijabs*g_rijabs + &
      d_potential_epp(epp(ids, jds, idj), potentialf(ui, uj, rij))*d_epp_anisotropic*(&
      d_anisotropic_ids(ui, uj, urij, chiepsilon)*g_ids(ui, urij, rijabs) + &
      d_anisotropic_ids(uj, ui, urij, chiepsilon)*g_ids(uj, urij, rijabs)) + &
      d_potential_gb_R(gb_epsilon(ui, uj, urij), gb_R(ui, uj, rij))*d_gb_R_sigma()*&
      d_sigma_anisotropic(anisotropic(ui, uj, urij, chisigma))*( &
      d_anisotropic_ids(ui, uj, urij, chisigma)*g_ids(ui, urij, rijabs) + &
      d_anisotropic_ids(uj, ui, urij, chisigma)*g_ids(uj, urij, rijabs))
  end function


  !! Calculates the interaction energy of a GB particle (i) and a Lennard-Jones
  !! particle (j)
  pure subroutine gblj_potential(ui, rij, energy, overlap)
    real(dp), intent(in) :: ui(3), rij(3)
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    real(dp) :: rijabs, urij(3), sigma, epsilon, R
    rijabs = sqrt(dot_product(rij, rij))
    urij = rij/rijabs
    sigma = gblj_sigma_0 / sqrt(1 - chisigma * dot_product(ui, urij)**2)
    epsilon = gblj_epsilon_0 * (1 - chiepsilon * dot_product(ui, urij)**2)
    R = (rijabs - sigma + gblj_sigma_0) / gblj_sigma_0
    if (hardcore > R) then 
      overlap = .true.
      energy = 0._dp
    else
      overlap = .false.
      energy = 4 * epsilon * R**(-6) * (R**(-6) - 1)
    end if
  end subroutine



end module


