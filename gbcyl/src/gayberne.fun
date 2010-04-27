test_suite gayberne



test zero_at_cross_contact
  use nrtype, only: dp
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
  assert_real_equal(0._dp, e_gb)
  assert_false(overlap)
!  assert_real_equal(0._dp, potential(ui, uj, rij))
end test



test kappa_sigma_defined
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
  assert_real_equal(kappa_sigma, sigma(ui, uj, urij_ee)/sigma(ui, uj, urij_ss))
end test



test kappa_epsilon_defined
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
  assert_real_equal(kappa_epsilon, epsilon(ui, uj, urij_ss)/epsilon(ui, uj, urij_ee))  
end test



test reduces_to_lennard_jones
  use nrtype, only: dp
  real(dp) :: sigma_0 = 3.435
  real(dp) :: epsilon_0 = 234.234
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
  assert_real_equal(lj_potential, e_gb)
  assert_false(overlap)
end test



test small_separation
  use nrtype, only: sp, dp
  intrinsic huge
  intrinsic tiny
  real(dp) :: kappa_sigma = 4.4
  real(dp) :: kappa_epsilon = 20._dp
  real(dp) :: mu = 1._dp
  real(dp) :: nu = 1._dp
  real(dp) :: sigma_0 = 1._dp
  real(dp) :: epsilon_0 = 1._dp
  real(dp), dimension(3) :: rij 
  real(dp), dimension(3) :: ui = (/1.0, 0.0, 0.0/)
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
  assert_true(overlap)
end test

test normalop
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
  expected = 4._dp * (dot_product(rij, rij)**(-6) - dot_product(rij, rij)**(-3))
  khi = (kappa_sigma**2 - 1._dp)/(kappa_sigma**2 + 1._dp)
  expected = expected/sqrt(1._dp - khi**2)
  !write(*, *) 'expected = ', expected
  assert_real_equal(expected, e_gb)
end test

end test_suite