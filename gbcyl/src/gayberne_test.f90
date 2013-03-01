module gayberne_test
use gayberne
use ftnunit
use nrtype
implicit none
private 
public :: test_all

real(dp), parameter :: margin = 1.e-9_dp
!! The relative tolerance for reals to be comparable.
!! @see ftnunit documentation. Two values v1 and v2 are considered comparable 
!! if abs( v1 - v2 ) < margin * (abs(v1)+abs(v2)) / 2

real(dp), dimension(3), parameter :: ex = (/1._dp, 0._dp, 0._dp/), &
ey = (/0._dp, 1._dp, 0._dp/), ez = (/0._dp, 0._dp, 1._dp/)


contains

subroutine test_all
  call test(test_zero_at_cross_contact, "Calculate potential at cross contact")
  call test(test_kappa_sigma_defined, "Test that definition of kappa for the "//&
  "contact distance calculation is correct")
  call test(test_kappa_epsilon_defined, "Test that definition of kappa for "//&
  "the well-depth calculation is correct")
  call test(test_normalop, "Test normal operation.")
  call test(test_reduces_to_lennard_jones, "Test that the potential reduces "//&
  "to a Lennard-Jones potential.")
  call test(test_small_separation, "Test that small separation results as "//&
  "an overlap.")
  call test(test_ssderivative, "Test the ratio of potential and its "//&
  "derivative in side by side configuration")
  call test(test_zerosofforce, "Test the zeros of force in ss, ee and T configurations.")
  call test(test_forcevsfinitedifference, "Test the force against finite "//&
  "difference approximation")
end subroutine

!! 1. "end test" -> "end subroutine" 
!! 2. "test " -> "subroutine test_" 
!! 3.1. 'assert_real_equal(' -> 'assert_comparable('
!! 3.2.  ')' ', resolution, "text")'
!!
!!
subroutine test_zero_at_cross_contact
  real(dp) :: kappa_sigma = 4.4_dp
  real(dp) :: kappa_epsilon = 20._dp
  real(dp) :: mu = 1.235_dp
  real(dp) :: nu = 1.6546_dp
  real(dp) :: sigma_0 = 1.234_dp
  real(dp) :: epsilon_0 = 1.2798_dp
  real(dp), dimension(3) :: ui = (/1._dp, 0._dp, 0._dp/)
  real(dp), dimension(3) :: uj = (/0._dp, 1._dp, 0._dp/)
  real(dp), dimension(3) :: rij
  real(dp) :: e_gb
  logical :: overlap
  rij = (/0.0_dp, 0.0_dp, sigma_0/)
  call gayberne_init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  call potential(ui, uj, rij, e_gb, overlap)
  call assert_comparable(0._dp, real(e_gb, dp), margin, &
  'Zero at cross contact')
  call assert_false(overlap, 'Overlap is false')
!  call assert_comparable(0._dp, potential(ui, uj, rij))
end subroutine

subroutine test_kappa_sigma_defined
  use nrtype, only: dp
  real(dp) :: kappa_sigma = 4.4_dp
  real(dp) :: kappa_epsilon = 20._dp
  real(dp) :: mu = 5.235_dp
  real(dp) :: nu = 3.6546_dp
  real(dp) :: sigma_0 = 1.234_dp
  real(dp) :: epsilon_0 = 1.2798_dp
  real(dp), dimension(3) :: ui = (/1._dp, 0._dp, 0._dp/)
  real(dp), dimension(3) :: uj 
  real(dp), dimension(3) :: urij_ee
  real(dp), dimension(3) :: urij_ss = (/0._dp, 1._dp, 0._dp/) 
  uj = ui 
  urij_ee = ui
  call gayberne_init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  call assert_comparable(kappa_sigma, &
  sigma(ui, uj, urij_ee)/sigma(ui, uj, urij_ss), margin, &
  "kappasigma defined correctly")
end subroutine

subroutine test_kappa_epsilon_defined
  use nrtype, only: dp
  real(dp) :: kappa_sigma = 4.4_dp
  real(dp) :: kappa_epsilon = 20._dp
  real(dp) :: mu = 6.235_dp
  real(dp) :: nu = 23.6546_dp
  real(dp) :: sigma_0 = 1.234_dp
  real(dp) :: epsilon_0 = 1.2798_dp
  real(dp), dimension(3) :: ui = (/1._dp, 0._dp, 0._dp/)
  real(dp), dimension(3) :: uj 
  real(dp), dimension(3) :: urij_ee
  real(dp), dimension(3) :: urij_ss = (/0._dp, 1._dp, 0._dp/) 
  uj = ui 
  urij_ee = ui
  call gayberne_init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  call assert_comparable(kappa_epsilon, &
  gb_epsilon(ui, uj, urij_ss)/gb_epsilon(ui, uj, urij_ee), margin, &
  "kappaepsilon defined correctly.")  
end subroutine

subroutine test_reduces_to_lennard_jones
  use nrtype, only: dp
  real(dp) :: sigma_0 = 3.435_dp
  real(dp) :: epsilon_0 = 234.234_dp
  real(dp) :: lj_potential
  real(dp) :: lj_6
  real(dp), dimension(3) :: rij 
  real(dp), dimension(3) :: ui = (/1._dp, 0._dp, 0._dp/)
  real(dp), dimension(3) :: uj = (/0._dp, 1._dp, 0._dp/)
  real(dp) :: r_absolute
  real(dp) :: e_gb
  logical :: overlap
  r_absolute = sigma_0 + 1._dp
  rij = (/-r_absolute, 0._dp, 0._dp/)
  lj_6 = (sigma_0 / r_absolute)**6
  lj_potential = 4._dp * epsilon_0 * lj_6 * (lj_6 - 1._dp)
  call gayberne_init(1._dp, 1._dp, 6.234_dp, 2.84687_dp, sigma_0, epsilon_0)
  call potential(ui, uj, rij, e_gb, overlap)
  call assert_comparable(lj_potential, e_gb, margin, &
  "Reduces to Lennard-Jones.")
  call assert_false(overlap, "Overlap is false.")
end subroutine

subroutine test_ssderivative
  use nrtype, only: dp
  real(dp), parameter :: mu = 1.5_dp, nu=1.3_dp
  real(dp), parameter :: sigma0 = 1.2_dp, epsilon0 = 0.9_dp
  real(dp), parameter :: kappasigma = 4.4_dp, kappaepsilon=20._dp
  real(dp) :: deriv, energy
  logical :: overlap
  call gayberne_init(kappasigma, kappaepsilon, mu, nu, sigma0, epsilon0)
  deriv = d_potential(ez, ez, ex*sigma0, 1)
  call potential(ez, ez, ex*sigma0, energy, overlap)
  call assert_comparable(deriv/energy, -6._dp*sigma0, margin, "Side-by-side derivative ok.")
end subroutine

subroutine test_small_separation
  use nrtype, only: sp, dp
  intrinsic huge
  intrinsic tiny
  real(dp) :: kappa_sigma = 4.4_dp
  real(dp) :: kappa_epsilon = 20._dp
  real(dp) :: mu = 1._dp
  real(dp) :: nu = 1._dp
  real(dp) :: sigma_0 = 1._dp
  real(dp) :: epsilon_0 = 1._dp
  real(dp), dimension(3) :: rij 
  real(dp), dimension(3) :: ui = (/1._dp, 0._dp, 0._dp/)
  real(dp), dimension(3) :: uj
  real(dp) :: r_absolute
  real(dp) :: e_gb
  logical :: overlap
  real(dp) :: small_number
  real(dp) :: hard_core = 0.6_dp
  small_number = 1.e-9_dp
  r_absolute = hard_core - small_number
  call gayberne_init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  uj = ui
  rij = (/0.0_dp, r_absolute, 0.0_dp/)
  call potential(ui, uj, rij, e_gb, overlap)
  call assert_true(overlap, "Overlap at small separation.")
end subroutine

subroutine test_normalop
  use nrtype
  real(dp) :: kappa_sigma = 4.4_dp
  real(dp) :: kappa_epsilon = 20._dp
  real(dp) :: mu = 1._dp
  real(dp) :: nu = 1._dp
  real(dp) :: sigma_0 = 1._dp
  real(dp) :: epsilon_0 = 1._dp
  real(dp), dimension(3) :: rij = (/1._dp, 1._dp, 0._dp/)
  real(dp), dimension(3) :: ui = (/0._dp, 0._dp, 1._dp/)
  real(dp) :: r_absolute
  real(dp) :: e_gb
  real(dp) :: expected
  logical :: overlap
  real(dp) :: khi
  call gayberne_init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  call potential(ui, ui, rij, e_gb, overlap)
  expected = 4._dp * (dot_product(rij, rij)**(-6) - &
  dot_product(rij, rij)**(-3))
  khi = (kappa_sigma**2 - 1._dp)/(kappa_sigma**2 + 1._dp)
  expected = expected/sqrt(1._dp - khi**2)
  !write(*, *) 'expected = ', expected
  call assert_comparable(expected, e_gb, margin, &
  "Normal operation.")
end subroutine

subroutine test_zerosofforce
  use nrtype
  real(dp) :: kappa_sigma = 4.4_dp
  real(dp) :: kappa_epsilon = 20._dp
  real(dp) :: mu = 1._dp
  real(dp) :: nu = 1._dp
  real(dp) :: sigma_0 = 1._dp
  real(dp) :: epsilon_0 = 1._dp
  real(dp) :: r_absolute
  real(dp) :: e_gb
  real(dp) :: expected
  logical :: overlap
  real(dp) :: khi
  real(dp) :: force
  real(dp), parameter :: sep = 2._dp**(1._dp/6._dp)
  call gayberne_init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  !! side-by-side configuration
  force = d_potential(ez, ez, ex*sep*sigma_0, 1)
  call assert_true(abs(force) < margin, "Force zero at ss configuration, r=2**(1/8)sigma_0")
  !! end-to-end configuration
  force = d_potential(ez, ez, ez*sigma_0*(sep-1._dp+kappa_sigma), 3)
  call assert_true(abs(force) < margin, "Force zero at ee configuration, r=(2**(1/8)-1+kappa_sigma)*sigma_0")
  force = d_potential(ez, ex, ex*sigma_0*(sep-1._dp+sqrt((kappa_sigma**2+1._dp)/2._dp)),1)
  call assert_true(abs(force) < margin, "Force zero at T configuration, "//&
  "r=(2**(1/8)-1+sqrt((kappa_sigma**2+1)/2))*sigma_0")
end subroutine

subroutine test_forcevsfinitedifference
  use nrtype
  real(dp), parameter :: kappa_sigma = 4.4_dp
  real(dp), parameter :: kappa_epsilon = 20._dp
  real(dp), parameter :: mu = 1._dp
  real(dp), parameter :: nu = 1._dp
  real(dp), parameter :: sigma_0 = 1._dp
  real(dp), parameter :: epsilon_0 = 1._dp
  real(dp) :: r_absolute
  logical :: overlap
  real(dp), dimension(3), parameter :: r0 = ex*(kappa_sigma-1.0_dp)*2._dp/3._dp
  real(dp), parameter :: dr = 1e-6_dp
  real(dp), dimension(3), parameter :: r1 = r0 + ex*dr
  real(dp), dimension(3), parameter :: uj = (/sqrt(2._dp)**(-1), 0._dp, sqrt(2._dp)**(-1)/) 
  real(dp) :: force, energy1, energy0
  real(dp),dimension(3) :: gradient
  call gayberne_init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  call potential(ez, uj, r0, energy0, overlap)
  call potential(ez, uj, r1, energy1, overlap)
  force= d_potential(ez, uj, r0+0.5_dp*dr*ex, 1)
  call assert_comparable(force, (energy1-energy0)/dr, margin, "Force comparable "//&
  "to finite difference approximation.")
  gradient= g_potential(ez, uj, r0+0.5_dp*dr*ex)
  write(*, *) gradient
  call assert_comparable(gradient(1), &
  (energy1-energy0)/dr, margin, "New force comparable to finite difference "//&
  "approximation.")
end subroutine

end module 
