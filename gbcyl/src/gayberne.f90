module gayberne
  use nrtype, only: sp, dp
  use class_parameterizer
  use class_parameter_writer
  implicit none
  private

  public :: init
  public :: potential
  public :: sigma
  public :: sigma_0
  public :: kappa_sigma
  public :: epsilon
  public :: gb_write_parameters

  !! Parameters of the Gay-Berne potential. Documentation:
  !! @see Luckhurst & et.al J.Chem.Phys, Vol. 110, No. 14
  !!
  real(dp), save :: kappa_sigma_   = 4.4_dp   !! = sig_e/sig_s
  real(dp), save :: kappa_epsilon_ = 20._dp !! = eps_s/eps_e  
  real(dp), save :: mu_            = 1._dp
  real(dp), save :: nu_            = 1._dp
  real(dp), save :: sigma_0_       = 1._dp
  real(dp), save :: epsilon_0_     = 1._dp
  !! Parameters below are derived from parameters above.
  real(dp), save :: chi_epsilon_
  real(dp), save :: chi_sigma_
  real(dp), save :: chi_sigma_squared_

  interface potential
    module procedure potential_s
  end interface

  interface init
    module procedure init_parameterizer, init_old
  end interface

  contains

  subroutine init_parameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call get_parameter(reader, 'gb_kappa_sigma', kappa_sigma_)
    call get_parameter(reader, 'gb_kappa_epsilon', kappa_epsilon_)
    call get_parameter(reader, 'gb_mu', mu_)
    call get_parameter(reader, 'gb_nu', nu_)
    call get_parameter(reader, 'gb_sigma_0', sigma_0_)
    call get_parameter(reader, 'gb_epsilon_0', epsilon_0_)
    chi_epsilon_ = &
      (kappa_epsilon_**(1._dp / mu_) - 1._dp) / &
      (kappa_epsilon_**(1._dp / mu_)+ 1._dp)
    chi_sigma_ = (kappa_sigma_ * kappa_sigma_ - 1._dp) / &
      (kappa_sigma_ * kappa_sigma_ + 1._dp)
    chi_sigma_squared_ = chi_sigma_**2
  end subroutine

  subroutine init_old(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  implicit none
  real(dp), intent(in) :: kappa_sigma
  real(dp), intent(in) :: kappa_epsilon
  real(dp), intent(in) :: mu
  real(dp), intent(in) :: nu
  real(dp), intent(in) :: sigma_0
  real(dp), intent(in) :: epsilon_0
    kappa_sigma_ = kappa_sigma
    kappa_epsilon_ = kappa_epsilon
    mu_ = mu
    nu_ = nu
    sigma_0_ = sigma_0
    epsilon_0_ = epsilon_0
    chi_epsilon_ = &
      (kappa_epsilon_**(1._dp / mu_) - 1._dp) / &
      (kappa_epsilon_**(1._dp / mu_)+ 1._dp)
    chi_sigma_ = (kappa_sigma_ * kappa_sigma_ - 1._dp) / &
      (kappa_sigma_ * kappa_sigma_ + 1._dp)
    chi_sigma_squared_ = chi_sigma_**2
  end subroutine

  subroutine gb_write_parameters(writer)
    type(parameter_writer), intent(in) :: writer
    call write_comment(writer, 'Gay-Berne potential parameters')
    call write_parameter(writer, 'gb_kappa_sigma', kappa_sigma_)
    call write_parameter(writer, 'gb_kappa_epsilon', kappa_epsilon_)
    call write_parameter(writer, 'gb_mu', mu_)
    call write_parameter(writer, 'gb_nu', nu_)
    call write_parameter(writer, 'gb_sigma_0', sigma_0_)
    call write_parameter(writer, 'gb_epsilon_0', epsilon_0_)    
  end subroutine

  !! Calculates the Gay-Berne potential for two particles 
  !! particlei and particlej. If r_gb=rij-sig+sig0<=0.6 
  !! routine returns with ovrlp=.true. 
  !!
  !! :TODO: check if ovrlp is really needed and in which conditions.
  !! :TODO: make hard_core adjustable
  !!
  pure subroutine potential_s(ui, uj, rij, gbV, ovrlp)
    implicit none
    real(dp), dimension(3), intent(in) :: rij, ui, uj
    real(dp), intent(out) :: gbV
    logical, intent(out) :: ovrlp
    real(dp) :: r_gb, gb6
    real(dp), parameter :: hard_core = 0.6_dp
    real(dp), dimension(3) :: urij
    gbV = 0._dp
    ovrlp = .false.
    r_gb = separation(ui, uj, rij)
    if(r_gb < hard_core) then
      ovrlp = .true.
    else
      gb6 = r_gb**(-6)
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
    separation = (r - sigma(ui, uj, urij) + sigma_0_) / sigma_0_
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
    sigma = sigma_help(ids, jds, idj)
  end function sigma

  !! :TODO: change this to take unit vectors as arguments. 
  real(dp) pure function sigma_help(ids, jds, idj)
  implicit none
  intrinsic sqrt
  real(dp), intent(in) :: ids, jds, idj
  real(dp) :: idssq, jdssq, idjsq
    idssq = ids*ids
    jdssq = jds*jds
    idjsq = idj*idj
    sigma_help = 1._dp - chi_sigma_ * (idssq + jdssq - 2._dp * chi_sigma_ * &
      ids * jds * idj) / (1._dp - chi_sigma_squared_ * idjsq)
    sigma_help = sigma_0_ / sqrt(sigma_help)
  end function sigma_help

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
    epsilon = epsilon_0_ * ep(idj)**nu_ * epp(ids, jds, idj)**mu_
  end function epsilon

  !! :TODO: change this to take unit vectors as arguments. 
  real(dp) pure function ep(idj)
  implicit none
  intrinsic sqrt
  real(dp), intent(in) :: idj
    ep = 1._dp / sqrt(1._dp - chi_sigma_squared_ * idj**2)
  end function ep

  !! :TODO: change this to take unit vectors as arguments.   
  real(dp) pure function epp(ids, jds, idj)
  implicit none
  real(dp),intent(in) :: ids, jds, idj
    epp = 1._dp - chi_epsilon_ * (ids**2 + jds**2 - 2._dp * chi_epsilon_ * &
      ids * jds * idj) / (1._dp - chi_epsilon_**2 * idj**2)
  end function epp

  function sigma_0()
  real(dp) :: sigma_0
    sigma_0 = sigma_0_
  end function sigma_0

  function kappa_sigma()
  real(dp) :: kappa_sigma
    kappa_sigma = kappa_sigma_
  end function kappa_sigma    

end module gayberne
