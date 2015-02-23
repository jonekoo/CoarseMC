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
module gblj
use nrtype
use class_parameterizer
use class_parameter_writer
implicit none

public :: gblj_init
public :: gblj_writeparameters
public :: gblj_potential
public :: gblj_force
public :: gblj_get_sigma_0
public :: gblj_r
public :: test_defaults

private

!> Well-depth when the LJ particle is on the side of the GB particle.
real(dp), save :: epsilon_0 = 1.55

!> The distance at which the GB-LJ potential crosses zero when the LJ
!! particle is on the side of the GB particle.
real(dp), save :: sigma_0 = 0.8 !! 3.6_dp / sigma0_aengstroms

!> Used for fine-tuning the relation of the potential energies for
!! different LJ particle positions with respect to the GB particle.
real(dp), save :: mu = 0.35

!> The ratio of well-depths when the LJ particle is at the end (ee) or
!! on the side (es) of the GB particle:
real(dp), save :: ee_to_es = 0.16
   
!> The ratio of contact distances when LJ is at the end (se) or on the
!! side (ss) of the GB particle:
real(dp), save :: se_to_ss = 3.32

real(dp), save :: chisigma
real(dp), save :: chiepsilon

!> Sets the GB-LJ distance at which the particles overlap.
real(dp), save :: hardcore = 0.3

!> Initializes the module.
interface gblj_init
   module procedure gblj_readparameters, gblj_init
end interface

contains


!> Calculates and saves chisigma and chiepsilon for later use.
subroutine gblj_init_internal()
  chisigma = 1 - se_to_ss ** (-2)
  chiepsilon = 1 - ee_to_es ** (1 / mu)
end subroutine


!> Initializes the module using a @p reader object, which is
!! responsible for reading the parameters e.g. from a file.
!! Parameters not given by @p reader are replaced by previously
!! set / default values.
subroutine gblj_readparameters(reader)
  type(parameterizer), intent(in) :: reader
  logical :: is_self_test = .false.
  call getparameter(reader, 'is_self_test', is_self_test)
  if (is_self_test) call test_defaults()
  call getparameter(reader, 'gblj_epsilon_0', epsilon_0)
  call getparameter(reader, 'gblj_sigma_0', sigma_0)
  call getparameter(reader, 'gblj_ee_to_es', ee_to_es)
  call getparameter(reader, 'gblj_se_to_ss', se_to_ss)
  call getparameter(reader, 'gblj_mu', mu)
  call getparameter(reader, 'gblj_hardcore', hardcore)
  call gblj_init_internal()
end subroutine


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
subroutine gblj_init(gblj_epsilon_0, gblj_sigma_0, gblj_ee_to_es, &
     gblj_se_to_ss, gblj_mu, gblj_hardcore)
  real(dp), intent(in), optional :: gblj_epsilon_0
  real(dp), intent(in), optional :: gblj_sigma_0
  real(dp), intent(in), optional :: gblj_ee_to_es
  real(dp), intent(in), optional :: gblj_se_to_ss
  real(dp), intent(in), optional :: gblj_mu
  real(dp), intent(in), optional :: gblj_hardcore
  if (present(gblj_epsilon_0)) epsilon_0 = gblj_epsilon_0
  if (present(gblj_sigma_0))   sigma_0 = gblj_sigma_0
  if (present(gblj_ee_to_es))  ee_to_es = gblj_ee_to_es
  if (present(gblj_se_to_ss))  se_to_ss = gblj_se_to_ss
  if (present(gblj_mu))        mu = gblj_mu
  if (present(gblj_hardcore))  hardcore = gblj_hardcore
  call gblj_init_internal()
end subroutine


!> Writes the parameters of this module to a file using the format and
!! output unit defined by @p writer.
subroutine gblj_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call writecomment(writer, &
       'Gay-Berne - Lennard-Jones potential parameters.')
  call writeparameter(writer, 'gblj_epsilon_0', epsilon_0)
  call writeparameter(writer, 'gblj_sigma_0', sigma_0)
  call writeparameter(writer, 'gblj_mu', mu)
  call writeparameter(writer, 'gblj_ee_to_es', ee_to_es)
  call writeparameter(writer, 'gblj_se_to_ss', se_to_ss)
  call writeparameter(writer, 'gblj_hardcore', hardcore)
end subroutine


!> Calculates the interaction energy of a GB particle (i) and a Lennard-Jones
!! particle (j)
!! 
!! @param ui      unit vector of orientation for the GB particle.
!! @param rij     vector between the centers of the two particles.
!! @param energy  interaction energy.
!! @param overlap indicates an overlap of particles when the potential
!!        would have a very high value.
!!
pure subroutine gblj_potential(ui, rij, energy, overlap)
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  real(dp) :: rijabs, urij(3), r
  rijabs = sqrt(dot_product(rij, rij))
  urij = rij/rijabs
  r = gblj_r(ui, rij)
  if (hardcore > r) then 
    overlap = .true.
    energy = 0._dp
  else
    overlap = .false.
    energy = 4 * gblj_epsilon(urij, ui) * (sigma_0 / r) ** 6 * &
         ((sigma_0 / r) ** 6 - 1)
  end if
end subroutine


!> Returns the force acting on the LJ particle (j) by the GB particle i.
pure function gblj_force(ui, rij) result(f)
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: f(3)
  real(dp) :: urij(3)
  urij = rij / sqrt(dot_product(rij, rij))
  f = 4. * sigma_0**6 / gblj_r(ui, rij)**13 * &
       ((gblj_r(ui, rij)**7 - gblj_r(ui, rij) * sigma_0) * &
       gblj_grad_eps(ui, rij) - 6. * gblj_epsilon(urij, ui) * &
       (gblj_r(ui, rij)**6 - 2. * sigma_0**6) * gblj_grad_r(ui, rij))
end function gblj_force

!> Returns the gradient of the well-depth gblj_epsilon. @p ui is the unit
!! orientation vector of the GB particle. @p rij is the vector from the
!! GB particle to the LJ particle.
pure function gblj_grad_eps(ui, rij) result(g)
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: g(3)
  g = -chiepsilon * (1 - gblj_aij(ui, rij) * chiepsilon)**(mu - 1) * mu * &
       epsilon_0 * gblj_grad_aij(ui, rij)
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

!> Returns the gradient of the distance function gblj_r. @p ui is the
!! unit vector of orientation of the GB particle. @p rij is the vector
!! from the GB particle to the LJ particle.
pure function gblj_grad_r(ui, rij) result(g)
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: g(3)
  real(dp) :: urij(3)
  urij = rij / sqrt(dot_product(rij, rij))
  g = urij - gblj_grad_sigma(ui, rij)
end function gblj_grad_r

!> Returns the gradient of the range function gblj_sigma. @p ui is the
!! unit orientation vector of the GB particle. @p rij is the vector
!! from the GB particle to the LJ particle. 
pure function gblj_grad_sigma(ui, rij) result(g)
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: g(3)
  g = chisigma * sigma_0 / (2. * sqrt(1. - gblj_aij(ui, rij) * chisigma)**3) &
       * gblj_grad_aij(ui, rij)
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
pure function gblj_sigma(urij, ui)
  real(dp), intent(in) :: urij(3), ui(3)
  real(dp) :: gblj_sigma
  gblj_sigma = sigma_0 / sqrt(1. - chisigma * dot_product(urij, ui)**2)
end function

!> Returns the anisotropic well-depth of the GB-LJ potential.
!! 
!! @param urij the unit vector pointing from the center of the GB
!!        particle to the center of the LJ particle.
!! @param ui the unit orientation vector of the GB particle. 
!!
pure function gblj_epsilon(urij, ui)
  real(dp), intent(in) :: urij(3), ui(3)
  real(dp) :: gblj_epsilon
  gblj_epsilon = epsilon_0 * (1. - chiepsilon * dot_product(urij, ui)**2)**mu
end function


!> The anisotropic distance function in the GB-LJ potential.
!! 
!! @param ui the unit orientation vector of the GB particle.
!! @param rij the vector from the center of the GB particle to the
!!        center of the LJ particle.
!!
pure function gblj_r(ui, rij)
  real(dp), intent(in) :: ui(3)
  real(dp), intent(in) :: rij(3)
  real(dp) :: gblj_r
  real(dp) :: urij(3)
  urij = rij / sqrt(dot_product(rij, rij))
  gblj_r = sqrt(dot_product(rij, rij)) - gblj_sigma(urij, ui) + sigma_0 
end function


!> Returns the contact distance of the GB-LJ potential when the LJ
!! particle is on the side of the GB particle.
pure function gblj_get_sigma_0()
  real(dp) :: gblj_get_sigma_0
  gblj_get_sigma_0 = sigma_0
end function

subroutine test_defaults()
  use m_fileunit
  real(dp) :: rij(3), ui(3), u
  integer :: dataunit
  real(dp) :: energy
  logical :: overlap
  integer :: ios
  call gblj_init_internal()
  dataunit = fileunit_getfreeunit()
  open(file = "gblj_testData.dat", unit = dataunit, action = 'READ')
  do 
    read(dataunit, fmt=*, iostat=ios) ui, rij, u
    if (ios /= 0) exit
    call gblj_potential(ui, rij, energy, overlap)
    if (.not. overlap .and. abs(u-energy) > 1.e-6_dp) then
      write(*, *) "gblj:test_defaults: ui=", ui, "rij=", rij, "abs(u-energy)=", abs(u-energy)
      stop
    end if
  end do
end subroutine

end module
