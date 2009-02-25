test_suite gayberne



test zero_at_cross_contact
  use nrtype, only: dp
  real(dp) :: kappa_sigma = 4.4
  real(dp) :: kappa_epsilon = 20.0
  real(dp) :: mu = 1.235
  real(dp) :: nu = 1.6546
  real(dp) :: sigma_0 = 1.234
  real(dp) :: epsilon_0 = 1.2798
  real(dp), dimension(3) :: ui = (/1.0_dp, 0.0_dp, 0.0_dp/)
  real(dp), dimension(3) :: uj = (/0.0_dp, 1.0_dp, 0.0_dp/)
  real(dp), dimension(3) :: rij
  rij = (/0.0_dp, 0.0_dp, sigma_0/)
  call init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  assert_real_equal(0.0_dp, potential(ui, uj, rij))
end test



test kappa_sigma_defined
  use nrtype, only: dp
  real(dp) :: kappa_sigma = 4.4
  real(dp) :: kappa_epsilon = 20.0
  real(dp) :: mu = 5.235
  real(dp) :: nu = 3.6546
  real(dp) :: sigma_0 = 1.234
  real(dp) :: epsilon_0 = 1.2798
  real(dp), dimension(3) :: ui = (/1.0_dp, 0.0_dp, 0.0_dp/)
  real(dp), dimension(3) :: uj 
  real(dp), dimension(3) :: urij_ee
  real(dp), dimension(3) :: urij_ss = (/0.0_dp, 1.0_dp, 0.0_dp/) 
  uj = ui 
  urij_ee = ui
  call init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  assert_real_equal(kappa_sigma, sigma(ui, uj, urij_ee)/sigma(ui, uj, urij_ss))
end test



test kappa_epsilon_defined
  use nrtype, only: dp
  real(dp) :: kappa_sigma = 4.4
  real(dp) :: kappa_epsilon = 20.0
  real(dp) :: mu = 6.235
  real(dp) :: nu = 23.6546
  real(dp) :: sigma_0 = 1.234
  real(dp) :: epsilon_0 = 1.2798
  real(dp), dimension(3) :: ui = (/1.0_dp, 0.0_dp, 0.0_dp/)
  real(dp), dimension(3) :: uj 
  real(dp), dimension(3) :: urij_ee
  real(dp), dimension(3) :: urij_ss = (/0.0_dp, 1.0_dp, 0.0_dp/) 
  uj = ui 
  urij_ee = ui
  call init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  assert_real_equal(kappa_epsilon, epsilon(ui, uj, urij_ss)/epsilon(ui, uj, urij_ee))  
end test



test reduces_to_lennard_jones
  use nrtype, only: dp
  real(dp) :: sigma_0 = 3.435
  real(dp) :: epsilon_0 = 234.234
  real(dp) :: lj_potential
  real(dp) :: lj_6
  real(dp), dimension(3) :: rij 
  real(dp), dimension(3) :: ui = (/1.0, 0.0, 0.0/)
  real(dp), dimension(3) :: uj = (/0.0, 1.0, 0.0/)
  real(dp) :: r_absolute
  r_absolute = sigma_0 + 1.0
  rij = (/-r_absolute, 0.0_dp, 0.0_dp/)
  lj_6 = (sigma_0/r_absolute)**(6)
  lj_potential = 4*epsilon_0*lj_6*(lj_6-1.0)
  call init(1.0_dp, 1.0_dp, 6.234_dp, 2.84687_dp, sigma_0, epsilon_0)
  assert_real_equal(lj_potential, potential(ui, uj, rij))
end test



test small_separation
  use nrtype, only: sp, dp
  intrinsic huge
  intrinsic tiny
  real(dp) :: kappa_sigma = 4.4
  real(dp) :: kappa_epsilon = 20.0
  real(dp) :: mu = 1.0
  real(dp) :: nu = 1.0
  real(dp) :: sigma_0 = 1.0
  real(dp) :: epsilon_0 = 1.0
  real(dp), dimension(3) :: rij 
  real(dp), dimension(3) :: ui = (/1.0, 0.0, 0.0/)
  real(dp), dimension(3) :: uj
  real(dp) :: r_absolute
  real(sp) :: single_precision_number
  r_absolute = 1e-6*tiny(single_precision_number)
  call init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  uj = ui
  rij = (/r_absolute, 0.0_dp, 0.0_dp/)
  assert_true(huge(r_absolute) > potential(ui, uj, rij))
  assert_true(sqrt(dot_product(rij, rij)) > tiny(r_absolute))
end test



end test_suite