module gayberne
  use nrtype, only: sp, dp
  use class_parameterizer
  use class_parameter_writer
  implicit none
  private

  public :: gayberne_init
  public :: potential
  public :: sigma
  public :: getsigma0
  public :: getkappasigma
  public :: epsilon
  public :: gb_writeparameters

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

  interface potential
    module procedure potentials
  end interface

  interface gayberne_init
    module procedure initparameterizer, initold
  end interface

  contains

  subroutine initparameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'gb_kappa_sigma', kappasigma)
    call getparameter(reader, 'gb_kappa_epsilon', kappaepsilon)
    call getparameter(reader, 'gb_mu', mu)
    call getparameter(reader, 'gb_nu', nu)
    call getparameter(reader, 'gb_sigma_0', sigma0)
    call getparameter(reader, 'gb_epsilon_0', epsilon0)
    chiepsilon = &
      (kappaepsilon**(1._dp / mu) - 1._dp) / &
      (kappaepsilon**(1._dp / mu)+ 1._dp)
    chisigma = (kappasigma * kappasigma - 1._dp) / &
      (kappasigma * kappasigma + 1._dp)
    chisigmasquared = chisigma**2
  end subroutine

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
    real(dp), parameter :: hardcore = 0.6_dp
    real(dp), dimension(3) :: urij
    gbV = 0._dp
    ovrlp = .false.
    rgb = separation(ui, uj, rij)
    if(rgb < hardcore) then
      ovrlp = .true.
    else
      gb6 = rgb**(-6)
      gbV = gb6 * (gb6 - 1._dp)
      urij = rij / sqrt(dot_product(rij, rij))
      gbV = 4._dp * epsilon(ui, uj, urij) * gbV
    end if
  end subroutine

  real(dp) pure function separation(ui, uj, rij)
    implicit none
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: rij
    real(dp) :: r
    real(dp), dimension(3) :: urij
    r = sqrt(dot_product(rij, rij))
    urij = rij/r
    separation = (r - sigma(ui, uj, urij) + sigma0) / sigma0
  end function separation

  real(dp) pure function sigma(ui, uj, urij)
    implicit none
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: urij
    real(dp) :: idj
    real(dp) :: ids
    real(dp) :: jds
    idj = dot_product(ui, uj) 
    ids = dot_product(ui, urij)
    jds = dot_product(uj, urij)
    sigma = sigmahelp(ids, jds, idj)
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

  real(dp) pure function epsilon(ui, uj, urij)
  implicit none
  real(dp), dimension(3), intent(in) :: ui
  real(dp), dimension(3), intent(in) :: uj
  real(dp), dimension(3), intent(in) :: urij
  real(dp) :: idj
  real(dp) :: ids
  real(dp) :: jds
    idj = dot_product(ui, uj) 
    ids = dot_product(ui, urij)
    jds = dot_product(uj, urij)
    epsilon = epsilon0 * ep(idj)**nu * epp(ids, jds, idj)**mu
  end function epsilon

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

end module
