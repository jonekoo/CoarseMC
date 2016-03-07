!> Module responsible for defining the computation of the (uniaxial)
!! Gay-Berne potential for the interaction of two rodlike molecules.
!!
!! @see J. G. Gay and B. J. Berne. J. Chem. Phys., 74(6):3316, 1981.
!! @see M. A. Bates and G. R. Luckhurst. J. Chem. Phys.,
!!      110(14):7087, 1999.
!!
module m_gayberne
  use num_kind
  use class_parameterizer
  use class_parameter_writer
  use json_module
  use m_json_wrapper, only: get_parameter
  implicit none

  !> Gay-Berne potential.
  !! 
  !! @see Luckhurst & et.al J.Chem.Phys, Vol. 110, No. 14
  !!
  type :: gayberne 
     !> Ratio of contact distances for end-to-end and side-by-side
     !! configurations of two particles.
     real(dp) :: kappasigma   = 4.4_dp   !! = sige/sigs
     
     !> Ratio of well-depths for side-by-side and end-to-end
     !! configurations of two particles.
     real(dp) :: kappaepsilon = 20._dp !! = epss/epse  
     
     !> This parameter is for fine-tuning how much the side-by-side
     !! configuration of two particles is favored as compared to the
     !! end-to-end configuration.  
     real(dp) :: mu            = 1._dp
     
     !> This parameter is for fine-tuning how much parallel alignment of
     !! two particles is favored as compared to other configurations
     !! (e.g. cross-configuration).
     real(dp) :: nu            = 1._dp
     
     !> The contact distance (where potential is zero) at
     !! cross-configuration of two particles.
     real(dp) :: sigma0       = 1._dp
     
     !> The well-depth at the cross-configuration of two particles.
     real(dp) :: epsilon0     = 1._dp
     
     
     !> Defines the GB distance gb_R, below which the particles overlap
     !! in gb_potential. Use to avoid numerical overflow.
     real(dp) :: hardcore = 0.6_dp
     
     !> Parameters below are derived from parameters above. These are
     !! saved to increase performance.
     real(dp) :: chiepsilon
     real(dp) :: chisigma
     real(dp) :: chisigmasquared
   
     !> Defines the procedure to be used for computing the well-depth
     !! function epsilon. This allows common special cases of the
     !! parameterization to be treated more efficiently than a single
     !! general function.
     procedure(gb_eps), pointer :: epsilon => null()
   contains
     procedure :: potential => gb_potential
     procedure :: force => gb_force
     procedure :: ep
     procedure :: epp
     procedure :: gb_R
     procedure :: potentialf
     procedure :: sigma
     procedure :: writeparameters => gb_writeparameters
     procedure :: to_json => gb_to_json
  end type gayberne

  !> Initializes the type
  interface gayberne
     module procedure initparameterizer, initold, gb_from_json
  end interface gayberne
  
  !> The interface for the well-depth function in the GB potential.
  interface
     pure function gb_eps(this, ui, uj, urij)
       import gayberne, dp
       implicit none
       class(gayberne), intent(in) :: this 
       real(dp), intent(in) :: ui(3), uj(3), urij(3)
       real(dp) :: gb_eps
     end function gb_eps
  end interface

contains

  !> Initializes the module using a parameterizer object.
  !! 
  !! @param[in] reader the parameterizer object which is responsible
  !! for getting the parameters for this module from some source, e.g.
  !! file. 
  !!
  function initparameterizer(reader) result(this)
    type(gayberne) :: this 
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'gb_kappa_sigma', this%kappasigma)
    call getparameter(reader, 'gb_kappa_epsilon', this%kappaepsilon)
    call getparameter(reader, 'gb_mu', this%mu)
    call getparameter(reader, 'gb_nu', this%nu)
    call getparameter(reader, 'gb_sigma_0', this%sigma0)
    call getparameter(reader, 'gb_epsilon_0', this%epsilon0)
    call getparameter(reader, 'gb_hardcore', this%hardcore)
    call init_common(this)
  end function initparameterizer

  !> Constructs a gayberne potential using the JSON in @p json_val.
  !! 
  function gb_from_json(json_val) result(this)
    type(json_value), pointer, intent(in) :: json_val
    type(gayberne) :: this 
    call get_parameter(json_val, 'gb_kappa_sigma', this%kappasigma, &
         error_lb=0._dp)
    call get_parameter(json_val, 'gb_kappa_epsilon', this%kappaepsilon, &
         error_lb=0._dp)
    call get_parameter(json_val, 'gb_mu', this%mu, error_lb=0._dp)
    call get_parameter(json_val, 'gb_nu', this%nu, error_lb=0._dp)
    call get_parameter(json_val, 'gb_sigma_0', this%sigma0, error_lb=0._dp)
    call get_parameter(json_val, 'gb_epsilon_0', this%epsilon0, error_lb=0._dp)
    call get_parameter(json_val, 'gb_hardcore', this%hardcore, error_lb=0._dp, &
         warn_ub=this%sigma0)
    call init_common(this)
  end function gb_from_json

  subroutine init_common(this)
    type(gayberne), intent(inout) :: this
    if (abs(this%mu - 1) < tiny(this%mu) .and. abs(this%nu - 1) < &
         tiny(this%nu)) then
       !write(*, *) '# gb_nu = 1 and gb_nu = 1. Using gb_epsilon_mu1nu1'
       this%epsilon => gb_epsilon_mu1nu1
    else
       this%epsilon => gb_epsilon_full
    end if
    this%chiepsilon = &
         (this%kappaepsilon**(1._dp / this%mu) - 1._dp) / &
         (this%kappaepsilon**(1._dp / this%mu)+ 1._dp)
    this%chisigma = (this%kappasigma * this%kappasigma - 1._dp) / &
         (this%kappasigma * this%kappasigma + 1._dp)
    this%chisigmasquared = this%chisigma**2
  end subroutine init_common


  !> Constructor.
  !! 
  !! @see M. A. Bates and G. R. Luckhurst, JCP 110(14), 7078, 1999 for
  !! a detailed discussion.
  !!
  !! @param kappasigma sets the axis ratio sigma_ee/sigma_ss of the
  !!        ellipsoidal molecule.
  !! @param kappaepsilon sets the ratio of well depths in
  !!        side-by-side and end-to-end configurations
  !!        epsilon_ss/epsilon_ee.
  !! @param mu adjusts how much side by side configuration is
  !!        favored.
  !! @param nu parameter for adjusting how much parallel alignment in
  !!        favored.
  !! @param sigma0 sets the contact distance for two ellipsoids in a
  !!        cross configuration. For a one-component Gay-Berne liquid
  !!        this can be set to 1.
  !! @param epsilon0 sets the well depth for the potential. For a one
  !!        component Gay-Berne liquid this can be set to 1.
  !! @param hardcore can be given to set a hard core to the potential.
  !!        In effect this defines the gb_R at which two particles
  !!        overlap.
  !! 
  function initold(kappasigma, kappaepsilon, mu, nu, sigma0, &
       epsilon0, hardcore) result(this)
    type(gayberne) :: this
    real(dp), intent(in) :: kappasigma
    real(dp), intent(in) :: kappaepsilon
    real(dp), intent(in) :: mu
    real(dp), intent(in) :: nu
    real(dp), intent(in) :: sigma0
    real(dp), intent(in) :: epsilon0
    real(dp), intent(in), optional :: hardcore
    this%kappasigma = kappasigma
    this%kappaepsilon = kappaepsilon
    this%mu = mu
    this%nu = nu
    this%sigma0 = sigma0
    this%epsilon0 = epsilon0
    if (present(hardcore)) this%hardcore = hardcore
    call init_common(this)
  end function initold


  !> Write the parameters of this module to the output unit and format
  !! defined by @p writer.
  subroutine gb_writeparameters(this, writer)
    class(gayberne), intent(in) :: this
    type(parameter_writer), intent(inout) :: writer
    call writecomment(writer, 'Gay-Berne potential parameters')
    call writeparameter(writer, 'gb_kappa_sigma', this%kappasigma)
    call writeparameter(writer, 'gb_kappa_epsilon', this%kappaepsilon)
    call writeparameter(writer, 'gb_mu', this%mu)
    call writeparameter(writer, 'gb_nu', this%nu)
    call writeparameter(writer, 'gb_sigma_0', this%sigma0)
    call writeparameter(writer, 'gb_epsilon_0', this%epsilon0)    
    call writeparameter(writer, 'gb_hardcore', this%hardcore)
  end subroutine gb_writeparameters


  !> Write the parameters of this module to the output unit and format
  !! defined by @p writer.
  subroutine gb_to_json(this, json_val)
    class(gayberne), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    call json_add(json_val, 'type', 'gayberne')
    call json_add(json_val, 'gb_kappa_sigma', this%kappasigma)
    call json_add(json_val, 'gb_kappa_epsilon', this%kappaepsilon)
    call json_add(json_val, 'gb_mu', this%mu)
    call json_add(json_val, 'gb_nu', this%nu)
    call json_add(json_val, 'gb_sigma_0', this%sigma0)
    call json_add(json_val, 'gb_epsilon_0', this%epsilon0)    
    call json_add(json_val, 'gb_hardcore', this%hardcore)
  end subroutine gb_to_json

  
  !> Calculates the Gay-Berne potential for two particles. If @p overlap
  !! is true, the potential is zero.
  !!
  !! @param ui,uj unit orientation vectors of the particles.
  !! @param rij vector from the center of particle i to particle j.
  !! @param energy the interaction energy.
  !! @param overlap true if the potential has a hard core and the particles
  !!        are too close to each other. 
  !!
  pure subroutine gb_potential(this, ui, uj, rij, energy, overlap)
    class(gayberne), intent(in) :: this
    real(dp), dimension(3), intent(in) :: rij, ui, uj
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    real(dp) :: rgb, gb6
    real(dp), dimension(3) :: urij
    energy = 0.0
    overlap = .false.
    if (all(rij == 0._dp)) then
       overlap = .true.
    else
       rgb = this%gb_R(ui, uj, rij)
       if(rgb < this%hardcore) then
          overlap = .true.
       else
          gb6 = rgb**(-6)
          energy = gb6 * (gb6 - 1._dp)
          urij = rij / sqrt(dot_product(rij, rij))
          energy = 4._dp * this%epsilon(ui, uj, urij) * energy
       end if
    end if
  end subroutine gb_potential
  

  !> Returns the Gay-Berne potential for two particles. Overlap of
  !! particles is not considered. 
  !!
  !! @param ui,uj the unit orientation vectors for the particles.
  !! @param rij the vector from the center of particle i to the center of
  !!        particle j.
  !! 
  real(dp) pure function potentialf(this, ui, uj, rij) result(energy)
    class(gayberne), intent(in) :: this
    real(dp), dimension(3), intent(in) :: rij, ui, uj
    real(dp) :: rgb, gb6
    real(dp), dimension(3) :: urij
    rgb = this%gb_R(ui, uj, rij)
    gb6 = rgb**(-6)
    energy = gb6 * (gb6 - 1._dp)
    urij = rij / sqrt(dot_product(rij, rij))
    energy = 4._dp * this%epsilon(ui, uj, urij) * energy
  end function potentialf


  !> Calculates the derivative of the potential with respect to the
  !! @p alpha coordinate of the vector @p rij.
  real(dp) pure function d_potential(this, ui, uj, rij, alpha)
    class(gayberne), intent(in) :: this
    real(dp), dimension(3), intent(in) :: rij, ui, uj
    integer, intent(in) :: alpha
    real(dp) :: rijabs
    real(dp), dimension(3) :: urij
    real(dp) :: energy
    logical :: overlap
    rijabs = sqrt(dot_product(rij, rij))
    urij = rij/rijabs
    call this%potential(ui, uj, rij, energy, overlap)
    d_potential = energy / this%epsilon(ui, uj, urij) * this%mu * & 
         rd_anisotropic2(ui, uj, urij, this%chiepsilon, alpha) / rijabs + &
         4._dp * this%epsilon(ui, uj, urij) * &
         (6._dp * this%gb_R(ui, uj, rij)**(-7) &
         - 12._dp * this%gb_R(ui, uj, rij)**(-13)) / &
         this%sigma0 * (rij(alpha) / rijabs + this%sigma0 / 2._dp / &
         sqrt(anisotropic(ui, uj, urij, this%chisigma))**3 * &
         rd_anisotropic2(ui, uj, urij, this%chisigma, alpha) / rijabs)
  end function d_potential
  
  real(dp) pure function rd_anisotropic2(ui, uj, urij, chi, alpha)
    real(dp), dimension(3), intent(in) :: ui, uj, urij
    real(dp), intent(in) :: chi
    integer, intent(in) :: alpha
    rd_anisotropic2 = 2 * urij(alpha) * &
         (1 - anisotropic(ui, uj, urij, chi)) - &
         chi / (1 - chi**2 * dot_product(ui,uj)**2) * &
         (tijalpha(ui, ui, urij, alpha) + tijalpha(uj, uj, urij, alpha) - &
         2 * chi * dot_product(ui,uj) * tijalpha(ui, uj, urij, alpha))
  end function rd_anisotropic2

  real(dp) pure function tijalpha(ui, uj, urij, alpha)
    real(dp), dimension(3), intent(in) :: ui, uj, urij
    integer, intent(in) :: alpha
    tijalpha = dot_product(uj, urij) * ui(alpha) + dot_product(ui, urij) * &
         uj(alpha)
  end function tijalpha
  
  !> The anisotropic distance function of the GB potential.
  !!
  !! @param ui,uj the unit orientation vectors of the two particles.
  !! @param rij the distance vector from the center of particle i to the center
  !!        of particle j.
  !!
  real(dp) pure function gb_R(this, ui, uj, rij)
    class(gayberne), intent(in) :: this
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: rij
    real(dp) :: r
    r = sqrt(dot_product(rij, rij))
    gb_R = (r - this%sigma(ui, uj, rij / r) + this%sigma0) / this%sigma0
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
  real(dp) pure function sigma(this, ui, uj, urij)
    class(gayberne), intent(in) :: this
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: urij
    sigma = this%sigma0 / sqrt(anisotropic(ui, uj, urij, this%chisigma))
  end function sigma

  real(dp) pure function sigmahelp(this, ids, jds, idj)
    class(gayberne), intent(in) :: this
    real(dp), intent(in) :: ids, jds, idj
    real(dp) :: idssq, jdssq, idjsq
    idssq = ids * ids
    jdssq = jds * jds
    idjsq = idj * idj
    sigmahelp = 1 - this%chisigma * (idssq + jdssq - &
         2 * this%chisigma * ids * jds * idj)&
         / (1 - this%chisigmasquared * idjsq)
    sigmahelp = this%sigma0 / sqrt(sigmahelp)
  end function sigmahelp

  real(dp) pure function gb_epsilon_full(this, ui, uj, urij)
    class(gayberne), intent(in) :: this
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: urij
    real(dp) :: idj
    idj = dot_product(ui, uj) 
    gb_epsilon_full = this%epsilon0 * this%ep(idj)**this%nu * &
         anisotropic(ui, uj, urij, this%chiepsilon)**this%mu
  end function gb_epsilon_full


  real(dp) pure function gb_epsilon_mu1nu1(this, ui, uj, urij)
    class(gayberne), intent(in) :: this
    real(dp), dimension(3), intent(in) :: ui
    real(dp), dimension(3), intent(in) :: uj
    real(dp), dimension(3), intent(in) :: urij
    real(dp) :: idj
    idj = dot_product(ui, uj) 
    gb_epsilon_mu1nu1 = this%epsilon0 * this%ep(idj) * &
         anisotropic(ui, uj, urij, this%chiepsilon)
  end function gb_epsilon_mu1nu1


  real(dp) pure function ep(this, idj)
    class(gayberne), intent(in) :: this
    real(dp), intent(in) :: idj
    ep = 1._dp / sqrt(1._dp - this%chisigmasquared * idj**2)
  end function ep

  real(dp) pure function epp(this, ids, jds, idj)
    class(gayberne), intent(in) :: this
    real(dp), intent(in) :: ids, jds, idj
    epp = 1 - this%chiepsilon * (ids**2 + jds**2 - &
         2 * this%chiepsilon * ids * jds * idj)&
         / (1 - this%chiepsilon**2 * idj**2)
  end function epp
  
  !> Returns the contact distance for two GB particles in a cross-configuration.
  real(dp) function getsigma0(this)
    class(gayberne), intent(in) :: this
    getsigma0 = this%sigma0
  end function getsigma0
  
  !> Returns the ratio of contact distances in the end-to-end and
  !! side-by-side configurations. 
  real(dp) function getkappasigma(this)
    class(gayberne), intent(in) :: this
    getkappasigma = this%kappasigma
  end function getkappasigma
  
  !! below here the "new" functions for force calculation
  
  real(dp) pure function d_anisotropic_ids(ui, uj, urij, chi)
    real(dp), dimension(3), intent(in) :: ui, uj, urij
    real(dp), intent(in) :: chi
    d_anisotropic_ids = -2._dp*chi/(1._dp-chi**2*dot_product(ui, uj)**2)*&
         (dot_product(ui, urij)-chi*dot_product(uj, urij)*dot_product(ui, uj))
  end function d_anisotropic_ids

  real(dp) pure function d_sigma_anisotropic(this, aniso)
    class(gayberne), intent(in) :: this
    real(dp), intent(in) :: aniso
    d_sigma_anisotropic = -0.5_dp * this%sigma0 * sqrt(aniso)**(-3)
  end function d_sigma_anisotropic
  
  real(dp) pure function d_gb_R_sigma(sigma0)
    real(dp), intent(in) :: sigma0
    d_gb_R_sigma = -1._dp/sigma0
  end function d_gb_R_sigma
  
  real(dp) pure function d_potential_gb_R(gb_eps, gb_r)
    real(dp), intent(in) :: gb_eps
    real(dp), intent(in) :: gb_r
    d_potential_gb_R = 4._dp*gb_eps*(6._dp*gb_r**(-7)-12._dp*gb_r**(-13))
  end function d_potential_gb_R

  real(dp) pure function d_potential_epp(epp, pot, mu)
    real(dp), intent(in) :: epp
    real(dp), intent(in) :: pot
    real(dp), intent(in) :: mu    
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
  pure function g_potential(this, ui, uj, rij)
    class(gayberne), intent(in) :: this
    real(dp), dimension(3), intent(in) :: ui, uj, rij
    real(dp) :: rijabs
    real(dp), dimension(3) :: g_potential
    real(dp) :: d_gb_R_rijabs 
    real(dp), dimension(3) :: urij
    real(dp) :: ids, jds, idj
    real(dp), dimension(3) :: g_rijabs
    real(dp), parameter :: d_epp_anisotropic = 1._dp 
    d_gb_R_rijabs = 1.0 / this%sigma0
    rijabs = sqrt(dot_product(rij, rij))
    urij = rij / rijabs
    g_rijabs = urij
    ids = dot_product(ui, urij)
    jds = dot_product(uj, urij)
    idj = dot_product(ui, uj)
    g_potential = d_potential_gb_R(this%epsilon(ui, uj, urij), &
         this%gb_R(ui, uj, rij)) &
         * d_gb_R_rijabs * g_rijabs + &
         d_potential_epp(this%epp(ids, jds, idj), &
         this%potentialf(ui, uj, rij), this%mu) * &
         d_epp_anisotropic * (d_anisotropic_ids(ui, uj, urij, this%chiepsilon)&
         * g_ids(ui, urij, rijabs) + &
         d_anisotropic_ids(uj, ui, urij, this%chiepsilon) &
         * g_ids(uj, urij, rijabs)) + &
         d_potential_gb_R(this%epsilon(ui, uj, urij), &
         this%gb_R(ui, uj, rij)) * &
         d_gb_R_sigma(this%sigma0) * &
         d_sigma_anisotropic(this, anisotropic(ui, uj, urij, this%chisigma)) * &
         ( &
         d_anisotropic_ids(ui, uj, urij, this%chisigma) * &
         g_ids(ui, urij, rijabs) + &
         d_anisotropic_ids(uj, ui, urij, this%chisigma) * &
         g_ids(uj, urij, rijabs))
  end function g_potential


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
  end function grad_anisotropic


  !> Returns the Gay-Berne force between two particles with unit
  !! orientation vectors @p ui,uj and center-to-center vector @p rij.
  pure function gb_force(this, ui, uj, rij)
    class(gayberne), intent(in) :: this
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
    e = this%epsilon(ui, uj, urij)
    grad_e = this%epsilon0 * this%ep(dot_product(ui, uj)) * &
         grad_anisotropic(ui, uj, rij, this%chiepsilon)
    
    !! The gradients of the Gay-Berne contact distance and distance
    grad_gb_sigma = -0.5 * this%sigma0 * &
         sqrt(anisotropic(ui, uj, urij, this%chisigma))**(-3) * &
         grad_anisotropic(ui, uj, rij, this%chisigma)
    grad_gb_R = (urij - grad_gb_sigma) / this%sigma0
    r = this%gb_R(ui, uj, rij)
    
    gb_force = -4 * (grad_e * (r**(-12) - r**(-6)) + &
         e * (6 * r**(-7) - 12 * r**(-13)) * grad_gb_R)
  end function gb_force
  
end module m_gayberne


