!> Implements calculation of the interaction energy for a 
!! Gay-Berne(GB)-Lennard-Jones(LJ) particle pair. 
!! 
!! @see J. Lintuvuori, M. Straka, and J. Vaara. Phys. Rev. E,
!!      75(3):031707, 2007.
!! @see H. Fukunaga, J. Takimoto, and M. Doi. J. Chem. Phys.,
!!      120(16):7792, 2004.
!! @see D. J. Cleaver, C. M. Care, M. P. Allen, and M. P. Neal.
!!      Phys. Rev. E, 54(1):559, 1996.
!! 
module m_gblj
  use num_kind
  use class_parameter_writer
  use json_module
  use m_json_wrapper, only: get_parameter
  implicit none
  
  public :: gblj_init
  public :: gblj_writeparameters
  public :: gblj_potential
  public :: gblj_force
  public :: gblj_get_sigma_0
  public :: gblj_r
  public :: test_defaults
  
  private
  
  
  type :: gblj_potential
     
     !> Well-depth when the LJ particle is on the side of the GB particle.
     real(dp) :: epsilon_0 = 1.55
     
     !> The distance at which the GB-LJ potential crosses zero when the LJ
     !! particle is on the side of the GB particle.
     real(dp) :: sigma_0 = 0.8 !! 3.6_dp / sigma0_aengstroms
     
     !> Used for fine-tuning the relation of the potential energies for
     !! different LJ particle positions with respect to the GB particle.
     real(dp) :: mu = 0.35
     
     !> The ratio of well-depths when the LJ particle is at the end (ee) or
     !! on the side (es) of the GB particle:
     real(dp) :: ee_to_es = 0.16
     
     !> The ratio of contact distances when LJ is at the end (se) or on the
     !! side (ss) of the GB particle:
     real(dp) :: se_to_ss = 3.32
     
     real(dp) :: chisigma
     real(dp) :: chiepsilon
     
     !> Sets the GB-LJ distance at which the particles overlap.
     real(dp) :: hardcore = 0.3
   contains
     
     procedure :: init_internal => gblj_init_internal
     procedure :: potential => gblj_value
     procedure :: sigma => gblj_sigma
     procedure :: epsilon => gblj_epsilon
     procedure :: force => gblj_force
     procedure :: grad_sigma => gblj_grad_sigma
     procedure :: grad_eps => gblj_grad_eps
     procedure :: r => gblj_r
     procedure :: grad_r => gblj_grad_r
     procedure :: writeparameters => gblj_writeparameters
     procedure :: to_json => gblj_to_json
  end type gblj_potential



  !> Initializes the module.
  interface gblj_potential
     module procedure gblj_init, gblj_from_json
  end interface gblj_potential
  
contains
  

!> Creates a gblj_potential from parameters in @p json_val.
function gblj_from_json(json_val) result(pot)
  type(json_value), pointer, intent(in) :: json_val
  type(gblj_potential) :: pot
  call get_parameter(json_val, 'gblj_epsilon_0', pot%epsilon_0, error_lb=0._dp)
  call get_parameter(json_val, 'gblj_sigma_0', pot%sigma_0, error_lb=0._dp)
  call get_parameter(json_val, 'gblj_ee_to_es', pot%ee_to_es, error_lb=0._dp)
  call get_parameter(json_val, 'gblj_se_to_ss', pot%se_to_ss, error_lb=0._dp)
  call get_parameter(json_val, 'gblj_mu', pot%mu, error_lb=0._dp)
  call get_parameter(json_val, 'gblj_hardcore', pot%hardcore, error_lb=0._dp, &
       warn_ub=pot%sigma_0)
  call pot%init_internal()
end function


!> Initializes this module. Parameters that are not given are replaced
!! with previously set/default values in the module. 
!!
!! @param gblj_epsilon_0 the well-depth for the GB-LJ potential.
!! @param gblj_sigma_0 the contact distance in a side-by-side
!!        configuration.
!! @param gblj_ee_to_es the ratio of well-depths in "end" and "side"
!!        configurations.
!! @param gblj_se_to_ss the ratio of contact distances in "end" and
!!        "side" configurations.
!! @param gblj_mu extra parameter for fine-tuning.
!! @param gblj_hardcore sets the gblj_r at which the GB and LJ particle
!!        overlap.
!!
function gblj_init(gblj_epsilon_0, gblj_sigma_0, gblj_ee_to_es, &
     gblj_se_to_ss, gblj_mu, gblj_hardcore) result(pot)
  real(dp), intent(in), optional :: gblj_epsilon_0
  real(dp), intent(in), optional :: gblj_sigma_0
  real(dp), intent(in), optional :: gblj_ee_to_es
  real(dp), intent(in), optional :: gblj_se_to_ss
  real(dp), intent(in), optional :: gblj_mu
  real(dp), intent(in), optional :: gblj_hardcore
  type(gblj_potential) :: pot
  if (present(gblj_epsilon_0)) pot%epsilon_0 = gblj_epsilon_0
  if (present(gblj_sigma_0))   pot%sigma_0 = gblj_sigma_0
  if (present(gblj_ee_to_es))  pot%ee_to_es = gblj_ee_to_es
  if (present(gblj_se_to_ss))  pot%se_to_ss = gblj_se_to_ss
  if (present(gblj_mu))        pot%mu = gblj_mu
  if (present(gblj_hardcore))  pot%hardcore = gblj_hardcore
  call pot%init_internal()
end function

!> Calculates and saves chisigma and chiepsilon for later use.
subroutine gblj_init_internal(this)
  class(gblj_potential), intent(inout) :: this
  this%chisigma = 1 - this%se_to_ss ** (-2)
  this%chiepsilon = 1 - this%ee_to_es ** (1 / this%mu)
end subroutine


!> Writes the parameters of this module to a file using the format and
!! output unit defined by @p writer.
subroutine gblj_writeparameters(this, writer)
  class(gblj_potential), intent(in) :: this
  class(parameter_writer), intent(inout) :: writer
  call writecomment(writer, &
       'Gay-Berne - Lennard-Jones potential parameters.')
  call writeparameter(writer, 'gblj_epsilon_0', this%epsilon_0)
  call writeparameter(writer, 'gblj_sigma_0', this%sigma_0)
  call writeparameter(writer, 'gblj_mu', this%mu)
  call writeparameter(writer, 'gblj_ee_to_es', this%ee_to_es)
  call writeparameter(writer, 'gblj_se_to_ss', this%se_to_ss)
  call writeparameter(writer, 'gblj_hardcore', this%hardcore)
end subroutine


subroutine gblj_to_json(this, json_val)
  class(gblj_potential), intent(in) :: this
  type(json_value), pointer, intent(inout) :: json_val
  call json_add(json_val, 'type', 'gblj')
  call json_add(json_val, 'gblj_epsilon_0', this%epsilon_0)
  call json_add(json_val, 'gblj_sigma_0', this%sigma_0)
  call json_add(json_val, 'gblj_mu', this%mu)
  call json_add(json_val, 'gblj_ee_to_es', this%ee_to_es)
  call json_add(json_val, 'gblj_se_to_ss', this%se_to_ss)
  call json_add(json_val, 'gblj_hardcore', this%hardcore)
end subroutine gblj_to_json


!> Calculates the interaction energy of a GB particle (i) and a Lennard-Jones
!! particle (j)
!! 
!! @param ui      unit vector of orientation for the GB particle.
!! @param rij     vector between the centers of the two particles.
!! @param energy  interaction energy.
!! @param overlap indicates an overlap of particles when the potential
!!        would have a very high value.
!!
pure subroutine gblj_value(this, ui, rij, energy, overlap)
  class(gblj_potential), intent(in) :: this
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  real(dp) :: rijabs, urij(3), r
  if (all(rij == 0)) then
     overlap = .true.
  else
     rijabs = sqrt(dot_product(rij, rij))
     urij = rij/rijabs
     r = this%r(ui, rij)
     if (this%hardcore > r) then 
        overlap = .true.
        energy = 0._dp
     else
        overlap = .false.
        energy = 4 * this%epsilon(urij, ui) * (this%sigma_0 / r) ** 6 * &
             ((this%sigma_0 / r) ** 6 - 1)
     end if
  end if
end subroutine gblj_value


!> Returns the force acting on the LJ particle (j) by the GB particle i.
pure function gblj_force(this, ui, rij) result(f)
  class(gblj_potential), intent(in) :: this
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: f(3)
  real(dp) :: urij(3)
  urij = rij / sqrt(dot_product(rij, rij))
  f = 4. * this%sigma_0**6 / this%r(ui, rij)**13 * &
       ((this%r(ui, rij)**7 - this%r(ui, rij) * this%sigma_0) * &
       this%grad_eps(ui, rij) - 6. * this%epsilon(urij, ui) * &
       (this%r(ui, rij)**6 - 2. * this%sigma_0**6) * this%grad_r(ui, rij))
end function gblj_force

!> Returns the gradient of the well-depth epsilon. @p ui is the unit
!! orientation vector of the GB particle. @p rij is the vector from the
!! GB particle to the LJ particle.
pure function gblj_grad_eps(this, ui, rij) result(g)
  class(gblj_potential), intent(in) :: this
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: g(3)
  g = -this%chiepsilon * (1 - gblj_aij(ui, rij) * this%chiepsilon)**(&
       this%mu - 1) * this%mu * this%epsilon_0 * gblj_grad_aij(ui, rij)
end function gblj_grad_eps

!> Returns the dot_product(ui, urij). @p ui is the unit orientation
!! vector of the GB particle. @p rij is the vector from the GB particle
!! to the LJ particle. urij is the unit vector with the same
!! orientation as @p rij
pure function gblj_aij(ui, rij) result(aij)
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: aij
  aij = dot_product(ui, rij) / sqrt(dot_product(rij, rij))
end function

!> Returns the gradient of the distance function r. @p ui is the
!! unit vector of orientation of the GB particle. @p rij is the vector
!! from the GB particle to the LJ particle.
pure function gblj_grad_r(this, ui, rij) result(g)
  class(gblj_potential), intent(in) :: this
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: g(3)
  real(dp) :: urij(3)
  urij = rij / sqrt(dot_product(rij, rij))
  g = urij - this%grad_sigma(ui, rij)
end function gblj_grad_r

!> Returns the gradient of the range function gblj_sigma. @p ui is the
!! unit orientation vector of the GB particle. @p rij is the vector
!! from the GB particle to the LJ particle. 
pure function gblj_grad_sigma(this, ui, rij) result(g)
  class(gblj_potential), intent(in) :: this
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: g(3)
  g = this%chisigma * this%sigma_0 / (2. * sqrt(1. - gblj_aij(ui, rij) * &
       this%chisigma)**3) * gblj_grad_aij(ui, rij)
end function gblj_grad_sigma

!> Returns the gradient of dot_product(ui, urij), where @p ui is the
!! unit orientation vector of the GB particle and urij is the unit
!! vector along @p rij, the vector from the GB particle to the LJ
!! particle. 
pure function gblj_grad_aij(ui, rij) result(g)
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: g(3)
  real(dp) :: urij(3)
  urij = rij / sqrt(dot_product(rij, rij))
  g = 2. * dot_product(ui, urij) * (-dot_product(ui, rij) / &
       sqrt(dot_product(rij, rij))**3 * rij + &
       ui / sqrt(dot_product(rij, rij))) 
end function gblj_grad_aij


!> Returns the anisotropic contact distance between a Gay-Berne and a
!! Lennard-Jones particle.
!!
!! @param urij the unit vector pointing from the center of the GB
!!        particle to the center of the LJ particle.
!! @param ui the unit vector of orientation for the GB particle.
!!
pure function gblj_sigma(this, urij, ui)
  class(gblj_potential), intent(in) :: this
  real(dp), intent(in) :: urij(3), ui(3)
  real(dp) :: gblj_sigma
  gblj_sigma = this%sigma_0 / sqrt(1. - this%chisigma * &
       dot_product(urij, ui)**2)
end function

!> Returns the anisotropic well-depth of the GB-LJ potential.
!! 
!! @param urij the unit vector pointing from the center of the GB
!!        particle to the center of the LJ particle.
!! @param ui the unit orientation vector of the GB particle. 
!!
pure function gblj_epsilon(this, urij, ui)
  class(gblj_potential), intent(in) :: this
  real(dp), intent(in) :: urij(3), ui(3)
  real(dp) :: gblj_epsilon
  gblj_epsilon = this%epsilon_0 * (1. - this%chiepsilon * &
       dot_product(urij, ui)**2)**this%mu
end function


!> The anisotropic distance function in the GB-LJ potential.
!! 
!! @param ui the unit orientation vector of the GB particle.
!! @param rij the vector from the center of the GB particle to the
!!        center of the LJ particle.
!!
pure function gblj_r(this, ui, rij)
  class(gblj_potential), intent(in) :: this
  real(dp), intent(in) :: ui(3)
  real(dp), intent(in) :: rij(3)
  real(dp) :: gblj_r
  real(dp) :: urij(3)
  urij = rij / sqrt(dot_product(rij, rij))
  gblj_r = sqrt(dot_product(rij, rij)) - this%sigma(urij, ui) + this%sigma_0 
end function


!> Returns the contact distance of the GB-LJ potential when the LJ
!! particle is on the side of the GB particle.
pure function gblj_get_sigma_0(this)
  class(gblj_potential), intent(in) :: this
  real(dp) :: gblj_get_sigma_0
  gblj_get_sigma_0 = this%sigma_0
end function

subroutine test_defaults(this)
  use m_fileunit
  class(gblj_potential), intent(inout) :: this
  real(dp) :: rij(3), ui(3), u
  integer :: dataunit
  real(dp) :: energy
  logical :: overlap
  integer :: ios
  call this%init_internal()
  dataunit = fileunit_getfreeunit()
  open(file = "gblj_testData.dat", unit = dataunit, action = 'READ')
  do 
    read(dataunit, fmt=*, iostat=ios) ui, rij, u
    if (ios /= 0) exit
    call this%potential(ui, rij, energy, overlap)
    if (.not. overlap .and. abs(u-energy) > 1.e-6_dp) then
       write(*, *) "gblj:test_defaults: ui=", ui, "rij=", rij, &
            "abs(u-energy)=", abs(u-energy)
      stop
    end if
  end do
end subroutine

end module m_gblj
