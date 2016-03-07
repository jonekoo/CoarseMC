!> Unit tests for the m_gayberne module.
module gayberne_pfunit
  use m_gayberne
  use pfunit
  use num_kind
  implicit none

  !> The relative tolerance for reals to be comparable.
  !! @see ftnunit documentation. Two values v1 and v2 are considered
  !! equal if abs( v1 - v2 ) < margin * (abs(v1)+abs(v2)) / 2
  real(dp), parameter :: margin = 1.e-9_dp
  
  real(dp), dimension(3), parameter :: ex = (/1._dp, 0._dp, 0._dp/), &
       ey = (/0._dp, 1._dp, 0._dp/), ez = (/0._dp, 0._dp, 1._dp/)
  
  
contains
  
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
    type(gayberne) :: gb
    rij = (/0.0_dp, 0.0_dp, sigma_0/)
    gb = gayberne(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
    call gb%potential(ui, uj, rij, e_gb, overlap)
    call assertEqual(0._dp, real(e_gb, dp), margin, &
         'Zero at cross contact')
    call AssertFalse(overlap, 'Overlap is false')
    !  call assertEqual(0._dp, potential(ui, uj, rij))
  end subroutine test_zero_at_cross_contact

  subroutine test_kappa_sigma_defined
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
    type(gayberne) :: gb
    uj = ui 
    urij_ee = ui
    gb = gayberne(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
    call assertEqual(kappa_sigma, &
         gb%sigma(ui, uj, urij_ee) / gb%sigma(ui, uj, urij_ss), margin, &
         "kappasigma defined correctly")
  end subroutine test_kappa_sigma_defined

  subroutine test_kappa_epsilon_defined
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
    type(gayberne) :: gb
    uj = ui 
    urij_ee = ui
    gb = gayberne(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
    call assertEqual(kappa_epsilon, &
         gb%epsilon(ui, uj, urij_ss) / gb%epsilon(ui, uj, urij_ee), &
         margin, "kappaepsilon defined correctly.")  
  end subroutine test_kappa_epsilon_defined

  subroutine test_reduces_to_lennard_jones
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
    type(gayberne) :: gb
    r_absolute = sigma_0 + 1._dp
    rij = (/-r_absolute, 0._dp, 0._dp/)
    lj_6 = (sigma_0 / r_absolute)**6
    lj_potential = 4._dp * epsilon_0 * lj_6 * (lj_6 - 1._dp)
    gb = gayberne(1._dp, 1._dp, 6.234_dp, 2.84687_dp, sigma_0, epsilon_0)
    call gb%potential(ui, uj, rij, e_gb, overlap)
    call assertEqual(lj_potential, e_gb, margin, &
         "Reduces to Lennard-Jones.")
    call AssertFalse(overlap, "Overlap is false.")
  end subroutine test_reduces_to_lennard_jones
  
  subroutine test_force_zeros
    use num_kind, only: dp
    real(dp), parameter :: mu = 1.5_dp, nu=1.3_dp
    real(dp), parameter :: sigma0 = 1.2_dp, epsilon0 = 0.9_dp
    real(dp), parameter :: kappasigma = 4.4_dp, kappaepsilon=20._dp
    real(dp) :: deriv(3), chisigma
    type(gayberne) :: gb
    gb = gayberne(kappasigma, kappaepsilon, mu, nu, sigma0, epsilon0)
    !! side-by-side:
    deriv = gb%force(ez, ez, ex * sigma0 * 2 ** (1._dp/6))
    call assertEqual(0._dp, deriv(1), margin, "Side-by-side derivative zero.")
    !! cross configuration:
    deriv = gb%force(ez, ey, ex * sigma0 * 2 ** (1._dp/6))
    call assertEqual(0._dp, deriv(1), margin, &
         "Cross configuration derivative zero.")
    !! T configuration:
    chisigma = (kappasigma * kappasigma - 1._dp) / &
         (kappasigma * kappasigma + 1._dp)
    deriv = gb%force(ez, ex, ex * sigma0 * &
         (2 ** (1._dp/6) - 1._dp/sqrt(1 - chisigma) + 1))
    call assertEqual(0._dp, deriv(1), margin, &
         "Cross configuration derivative zero.")
    !! End to end:
    deriv = gb%force(ez, ez, ez * sigma0 * &
         (2 ** (1._dp/6) - kappasigma + 1))
    call assertEqual(0._dp, deriv(1), margin, &
         "Cross configuration derivative zero.")
  end subroutine test_force_zeros
  
  subroutine test_small_separation
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
    type(gayberne) :: gb
    small_number = 1.e-9_dp
    r_absolute = hard_core - small_number
    gb = gayberne(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
    uj = ui
    rij = (/0.0_dp, r_absolute, 0.0_dp/)
    call gb%potential(ui, uj, rij, e_gb, overlap)
    call AssertTrue(overlap, "Overlap at small separation.")
  end subroutine test_small_separation

  subroutine test_normalop
    real(dp) :: kappa_sigma = 4.4_dp
    real(dp) :: kappa_epsilon = 20._dp
    real(dp) :: mu = 1._dp
    real(dp) :: nu = 1._dp
    real(dp) :: sigma_0 = 1._dp
    real(dp) :: epsilon_0 = 1._dp
    real(dp), dimension(3) :: rij = (/1._dp, 1._dp, 0._dp/)
    real(dp), dimension(3) :: ui = (/0._dp, 0._dp, 1._dp/)
    real(dp) :: e_gb
    real(dp) :: expected
    logical :: overlap
    real(dp) :: khi
    type(gayberne) :: gb
    gb = gayberne(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
    call gb%potential(ui, ui, rij, e_gb, overlap)
    expected = 4._dp * (dot_product(rij, rij)**(-6) - &
         dot_product(rij, rij)**(-3))
    khi = (kappa_sigma**2 - 1._dp)/(kappa_sigma**2 + 1._dp)
    expected = expected/sqrt(1._dp - khi**2)
    !write(*, *) 'expected = ', expected
    call assertEqual(expected, e_gb, margin, &
         "Normal operation.")
  end subroutine test_normalop
  
  subroutine test_forcevsfinitedifference
    real(dp), parameter :: kappa_sigma = 4.4
    real(dp), parameter :: kappa_epsilon = 20.
    real(dp), parameter :: mu = 1.
    real(dp), parameter :: nu = 1.
    real(dp), parameter :: sigma_0 = 1.
    real(dp), parameter :: epsilon_0 = 1.
    logical :: overlap
    real(dp), dimension(3), parameter :: r0 = ex * (kappa_sigma - 1.) * 2. / 3.
    real(dp), parameter :: dr = 1e-6
    !real(dp), dimension(3), parameter :: r1 = r0 + ex * dr
    real(dp), dimension(3) :: r1 
    real(dp), dimension(3), parameter :: uj = (/sqrt(2.)**(-1), 0., &
         sqrt(2.)**(-1)/) 
    real(dp) :: force(3), energy1, energy0
    real(dp),dimension(3) :: gradient
    type(gayberne) :: gb
    real(dp), parameter :: evecs(3, 3) = &
         reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
    integer :: i
    gb = gayberne(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)

    do i = 1, 3
       call gb%potential(ez, uj, r0, energy0, overlap)
       r1 = r0 + evecs(:, i) * dr
       call gb%potential(ez, uj, r1, energy1, overlap)
       force = gb%force(ez, uj, r0 + 0.5 * dr * evecs(:, i))
       call assertEqual(-(energy1 - energy0) / dr, force(i), margin, &
            "Force not comparable to finite difference approximation.")
       gradient= g_potential(gb, ez, uj, r0 + 0.5 * dr * evecs(:, i))
       call assertEqual(gradient(i), (energy1 - energy0) / dr, margin, &
            "New force comparable to finite difference approximation.")
    end do
end subroutine

end module 
