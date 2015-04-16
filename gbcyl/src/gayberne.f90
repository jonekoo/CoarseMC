module gayberne
  use nrtype, only: sp, dp
  implicit none

  public :: init
  public :: potential
  public :: sigma
  public :: epsilon
  public :: save_state
  public :: load_state

  private

  !! Parameters of the Gay-Berne potential. Documentation:
  !! @see Luckhurst & et.al J.Chem.Phys, Vol. 110, No. 14
  !!
  real(dp), save :: kappa_sigma_ !! = sig_e/sig_s
  real(dp), save :: kappa_epsilon_ !! = eps_s/eps_e  
  real(dp), save :: mu_
  real(dp), save :: nu_
  real(dp), save :: sigma_0_
  real(dp), save :: epsilon_0_

  real(dp), save :: chi_epsilon_
  real(dp), save :: chi_sigma_
  real(dp), save :: chi_sigma_squared_
  
  real(dp), save :: overlap_cutoff_

  namelist /gbgb_nml/ kappa_sigma_, kappa_epsilon_, mu_, nu_, sigma_0_, &
    epsilon_0_, chi_epsilon_, chi_sigma_, chi_sigma_squared_, &
    overlap_cutoff_
 
  interface potential
    module procedure potential
    module procedure potential_wt_overlap
  end interface


  contains



  subroutine init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
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
    overlap_cutoff_ = 0.6
    chi_epsilon_ = &
      & (kappa_epsilon_**(1.0/mu_) - 1.0)/(kappa_epsilon_**(1.0/mu_) + 1.0)
    chi_sigma_ = &
      & (kappa_sigma_*kappa_sigma_ - 1.0)/(kappa_sigma_*kappa_sigma_ + 1.0)
    chi_sigma_squared_ = chi_sigma_**2
  end subroutine init



  subroutine save_state(write_unit)
    implicit none
    integer, intent(in) :: write_unit
    write(write_unit, NML = gbgb_nml)
  end subroutine save_state



  subroutine load_state(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    read(read_unit, NML = gbgb_nml)
  end subroutine load_state



  !! Calculates the Gay-Berne potential for two particles 
  !! particlei and particlej. If r_gb=rij-sig+sig0<=0.6 
  !! routine returns with ovrlp=.true. 
  !!
  function potential_wt_overlap(ui, uj, rij, overlap) result(V)
    implicit none
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: rij
    logical, intent(out) :: overlap
    real(dp) :: V 
    if(separation(ui, uj, rij) < overlap_cutoff_) then
      overlap = .true.
      V = 0.0
    else
      V = potential(ui, uj, rij) 
    end if    
  end function potential_wt_overlap



  real(dp) function potential(ui, uj, rij)
    implicit none
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: rij
    real(dp) :: r_gb, gb6
    real(dp), dimension(3) :: urij
    r_gb = separation(ui, uj, rij)
    gb6 = r_gb**(-6)
    potential = gb6*(gb6-1.0)
    urij = rij/sqrt(dot_product(rij, rij))
    potential = 4*epsilon(ui, uj, urij)*potential
  end function



  real(dp) function separation(ui, uj, rij)
    implicit none
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: rij
    real(dp) :: r
    real(dp), dimension(3) :: urij
    r = sqrt(dot_product(rij, rij))
    urij = rij/r
    separation = (r - sigma(ui, uj, urij) + sigma_0_)/sigma_0_
  end function separation



  real(dp) function sigma(ui, uj, urij)
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
    sigma_help = 1-chi_sigma_*(idssq+jdssq-2*chi_sigma_*ids*jds*idj)/ &
      & (1-chi_sigma_squared_ * idjsq)
    sigma_help = sigma_0_/sqrt(sigma_help)
  end function sigma_help



  real(dp) function epsilon(ui, uj, urij)
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
    ep = 1.0/sqrt(1-chi_sigma_squared_*idj**2);
  end function ep



  !! :TODO: change this to take unit vectors as arguments.   
  real(dp) pure function epp(ids, jds, idj)
    implicit none
    real(dp),intent(in) :: ids, jds, idj
    epp = 1.0_dp-chi_epsilon_*(ids**2+jds**2-2.0_dp*chi_epsilon_*ids*jds*idj)/&
        (1.0_dp-chi_epsilon_**2 * idj**2)
  end function epp


      
end module gayberne
