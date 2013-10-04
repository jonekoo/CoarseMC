module gblj
use nrtype
implicit none

private

!! Parameters for the GB-LJ interaction:
real(dp), save :: epsilon_0 = 1._dp
real(dp), save :: sigma_0 = 1._dp
real(dp), save :: mu = 0.35

!! The ratio of well-depths when the LJ particle is at the end or on the side
!! of the GB particle:
real(dp), save :: ee_to_es = 0.16
   
!! The ratio of contact distances when LJ is at the end or on the side of the
!! GB particle:
real(dp), save :: se_to_ss = 3.32  

!! These are cached just for efficiency
real(dp), save :: chisigma
real(dp), save :: chiepsilon

contains

subroutine gblj_init(reader)
  call getparameter(reader, 'gblj_epsilon_0', epsilon_0)
  call getparameter(reader, 'gblj_sigma_0', sigma_0)
  call getparameter(reader, 'gblj_ee_to_es', ee_to_es)
  call getparameter(reader, 'gblj_se_to_ss', se_to_ss)
  call getparameter(reader, 'gblj_mu', mu)
  chisigma = 1._dp - se_to_ss ** (-2)
  chiepsilon = 1._dp - ee_to_es ** (1._dp / mu)
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
  real(dp) :: rijabs, urij(3), sigma, epsilon, r
  rijabs = sqrt(dot_product(rij, rij))
  urij = rij/rijabs
  r = gblj_r(ui, rij)
  if (hardcore > r) then 
    overlap = .true.
    energy = 0._dp
  else
    overlap = .false.
    energy = 4._dp * gblj_epsilon(urij, ui) * (gblj_sigma_0 / r) ** 6 * ((gblj_sigma_0 / r) ** 6 - 1._dp)
  end if
end subroutine

function gblj_force(ui, rij) 
  real(dp), intent(in) :: ui(3), rij(3)
  real(dp) :: gblj_force
  gblj_force =  g_gblj_epsilon(urij, ui) * (gblj_sigma_0 / r) ** 6 * &
   ((gblj_sigma_0 / r) ** 6 - 1._dp)
  gblj_force = gblj_force + gblj_epsilon(urij, ui) * (-6._dp) * &
    gblj_sigma_0 ** 6 / r ** 7 * (-2._dp * (gblj_sigma_0 / r) ** 6 + 1._dp) * &
    g_gblj_r(ui, rij)
  gblj_force = 4._dp * gblj_force
end function

!> Returns the anisotropic contact distance between a Gay-Berne and a Lennard-Jones particle.
!!
!! urij = unit vector between particles i and j.
!! ui   =  the unit vector of orientation for particle i.
!!
function gblj_sigma(urij, ui)
  real(dp), intent(in) :: urij(3), ui(3)
  real(dp) :: gblj_sigma
  gblj_sigma = sigma_0 / sqrt((1._dp - chisigma * dot_product(urij, ui) ** 2))
end function

function gblj_epsilon(urij, ui)
  real(dp), intent(in) :: urij(3), ui(3)
  real(dp) :: gblj_epsilon
  gblj_epsilon = epsilon_0 * (1._dp - chiepsilon * dot_product(urij, ui) ** 2) ** mu
end function

function gblj_r(ui, rij)
  real(dp), intent(in) :: ui(3)
  real(dp), intent(in) :: rij(3)
  real(dp) :: gblj_r
  real(dp) :: urij(3)
  urij = rij / sqrt(dot_product(rij, rij))
  gblj_r = sqrt(dot_product(rij, rij)) - gblj_sigma(urij, ui) + sigma_0 
end function

function get_gblj_sigma_0()
  real(dp) :: get_gblj_sigma_0
  get_gblj_sigma_0 = sigma_0
end function




end module
