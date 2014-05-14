module gblj
use nrtype
use class_parameterizer
use class_parameter_writer
use m_constants
implicit none

public :: gblj_init
public :: gblj_writeparameters
public :: gblj_potential
public :: gblj_force
public :: gblj_get_sigma_0
public :: gblj_r
public :: test_defaults

private

!! Parameters for the GB-LJ interaction:
real(dp), save :: epsilon_0 = 1.55_dp
real(dp), save :: sigma_0 = 3.6_dp / sigma0_aengstroms
real(dp), save :: mu = 0.35_dp

!! The ratio of well-depths when the LJ particle is at the end or on the side
!! of the GB particle:
real(dp), save :: ee_to_es = 0.16_dp
   
!! The ratio of contact distances when LJ is at the end or on the side of the
!! GB particle:
real(dp), save :: se_to_ss = 3.32_dp

!! These are cached just for efficiency
real(dp), save :: chisigma
real(dp), save :: chiepsilon

real(dp), save :: hardcore = 0.3_dp

contains

subroutine gblj_init_internal()
  chisigma = 1._dp - se_to_ss ** (-2)
  chiepsilon = 1._dp - ee_to_es ** (1._dp / mu)
end subroutine

subroutine gblj_init(reader)
  type(parameterizer), intent(in) :: reader
  logical :: is_self_test = .false.
  call getparameter(reader, 'is_self_test', is_self_test)
  if (is_self_test) call test_defaults()
  call getparameter(reader, 'gblj_epsilon_0', epsilon_0)
  call getparameter(reader, 'gblj_sigma_0', sigma_0)
  call getparameter(reader, 'gblj_ee_to_es', ee_to_es)
  call getparameter(reader, 'gblj_se_to_ss', se_to_ss)
  call getparameter(reader, 'gblj_mu', mu)
  call gblj_init_internal()
end subroutine

subroutine gblj_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call writeparameter(writer, 'gblj_epsilon_0', epsilon_0)
  call writeparameter(writer, 'gblj_sigma_0', sigma_0)
  call writeparameter(writer, 'gblj_mu', mu)
  call writeparameter(writer, 'gblj_ee_to_es', ee_to_es)
  call writeparameter(writer, 'gblj_se_to_ss', se_to_ss)
end subroutine

!> Calculates the interaction energy of a GB particle (i) and a Lennard-Jones
!! particle (j)
!! 
!! @p ui      = unit vector of orientation for the GB particle.
!! @p rij     = vector between the centers of the two particles.
!! @p energy  = interaction energy.
!! @p overlap = indicates an overlap of particles when the potential would have
!!              a very high value.
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
    energy = 4._dp * gblj_epsilon(urij, ui) * (sigma_0 / r) ** 6 * ((sigma_0 / r) ** 6 - 1._dp)
  end if
end subroutine

pure function gblj_force(ui, rij) 
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: gblj_force(3)
  real(dp) :: urij(3), r
  r = sqrt(dot_product(rij, rij))
  urij = rij / r
  gblj_force =  gblj_grad_epsilon(ui, urij, r) * (sigma_0 / gblj_r(ui, rij)) ** 6 * &
   ((sigma_0 / gblj_r(ui, rij)) ** 6 - 1._dp)
  gblj_force = gblj_force + gblj_epsilon(urij, ui) * (-6._dp) * &
    sigma_0 ** 6 / gblj_r(ui, rij) ** 7 * (-2._dp * (sigma_0 / gblj_r(ui, rij)) ** 6 + 1._dp) * &
    gblj_grad_r(ui, rij)
  gblj_force = -4._dp * gblj_force
end function

!> Returns the anisotropic contact distance between a Gay-Berne and a Lennard-Jones particle.
!!
!! urij = unit vector between particles i and j.
!! ui   =  the unit vector of orientation for particle i.
!!
pure function gblj_sigma(urij, ui)
  real(dp), intent(in) :: urij(3), ui(3)
  real(dp) :: gblj_sigma
  gblj_sigma = sigma_0 / sqrt((1._dp - chisigma * dot_product(urij, ui) ** 2))
end function

pure function gblj_epsilon(urij, ui)
  real(dp), intent(in) :: urij(3), ui(3)
  real(dp) :: gblj_epsilon
  gblj_epsilon = epsilon_0 * (1._dp - chiepsilon * dot_product(urij, ui) ** 2) ** mu
end function

pure function gblj_grad_es(se_0, ui, urij, rij, mu, chi)
  real(dp), intent(in) :: se_0, ui(3), urij(3), rij, mu, chi
  real(dp) :: gblj_grad_es(3)
  gblj_grad_es = -2._dp * se_0 * chi * mu * dot_product(ui, urij) * (urij * dot_product(ui, urij) + ui) / &
    (rij * (1._dp - chi * dot_product(ui, urij) ** 2) ** (1._dp - mu))  
end function

pure function gblj_grad_epsilon(ui, urij, rij)
  real(dp), intent(in) :: ui(3), urij(3), rij
  real(dp) :: gblj_grad_epsilon(3)
  gblj_grad_epsilon = gblj_grad_es(epsilon_0, ui, urij, rij, mu, chiepsilon)
end function

pure function gblj_grad_sigma(ui, urij, rij)
  real(dp), intent(in) :: ui(3), urij(3), rij
  real(dp) :: gblj_grad_sigma(3)
  gblj_grad_sigma = gblj_grad_es(sigma_0, ui, urij, rij, -0.5_dp, chisigma)
end function

pure function gblj_grad_r(ui, rij)
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: gblj_grad_r(3)
  real(dp) :: urij(3), r
  r = sqrt(dot_product(rij, rij))
  urij = rij / r
  gblj_grad_r = urij - gblj_grad_sigma(ui, urij, r)
end function 

pure function gblj_r(ui, rij)
  real(dp), intent(in) :: ui(3)
  real(dp), intent(in) :: rij(3)
  real(dp) :: gblj_r
  real(dp) :: urij(3)
  urij = rij / sqrt(dot_product(rij, rij))
  gblj_r = sqrt(dot_product(rij, rij)) - gblj_sigma(urij, ui) + sigma_0 
end function

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
