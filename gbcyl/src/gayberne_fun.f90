! gayberne_fun.f90 - a unit test suite for gayberne.f90
!
! funit generated this file from gayberne.fun
! at Thu Dec 04 15:59:26 +0200 2008

module gayberne_fun

 use gayberne

 implicit none

 logical :: noAssertFailed

 public :: test_gayberne

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0





 contains

 subroutine potential_matches_zero_at_cross_contact

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
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0.0_dp &
        +2*spacing(real(0.0_dp)) ) &
        .ge. &
        (potential(ui, uj, rij)) &
            .and. &
     (0.0_dp &
      -2*spacing(real(0.0_dp)) ) &
      .le. &
       (potential(ui, uj, rij)) )) then
      print *, " *Assert_Real_Equal failed* in test potential_matches_zero_at_cross_contact &
              &[gayberne.fun:18]"
      print *, "  ", "potential(ui, uj, rij) (", &
 potential(ui, uj, rij), &
  ") is not", &
 0.0_dp,&
 "within", &
  2*spacing(real(0.0_dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine potential_matches_zero_at_cross_contact




 subroutine kappa_sigma_matches_sigmaee_per_sigmass

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
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (kappa_sigma &
        +2*spacing(real(kappa_sigma)) ) &
        .ge. &
        (sigma(ui, uj, urij_ee)/sigma(ui, uj, urij_ss)) &
            .and. &
     (kappa_sigma &
      -2*spacing(real(kappa_sigma)) ) &
      .le. &
       (sigma(ui, uj, urij_ee)/sigma(ui, uj, urij_ss)) )) then
      print *, " *Assert_Real_Equal failed* in test kappa_sigma_matches_sigmaee_per_sigmass &
              &[gayberne.fun:38]"
      print *, "  ", "sigma(ui, uj, urij_ee)/sigma(ui, uj, urij_ss) (", &
 sigma(ui, uj, urij_ee)/sigma(ui, uj, urij_ss), &
  ") is not", &
 kappa_sigma,&
 "within", &
  2*spacing(real(kappa_sigma))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine kappa_sigma_matches_sigmaee_per_sigmass




 subroutine kappa_epsilon_matches_epsilonss_per_epsilonee

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
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (kappa_epsilon &
        +2*spacing(real(kappa_epsilon)) ) &
        .ge. &
        (epsilon(ui, uj, urij_ss)/epsilon(ui, uj, urij_ee)) &
            .and. &
     (kappa_epsilon &
      -2*spacing(real(kappa_epsilon)) ) &
      .le. &
       (epsilon(ui, uj, urij_ss)/epsilon(ui, uj, urij_ee)) )) then
      print *, " *Assert_Real_Equal failed* in test kappa_epsilon_matches_epsilonss_per_epsilonee &
              &[gayberne.fun:58]"
      print *, "  ", "epsilon(ui, uj, urij_ss)/epsilon(ui, uj, urij_ee) (", &
 epsilon(ui, uj, urij_ss)/epsilon(ui, uj, urij_ee), &
  ") is not", &
 kappa_epsilon,&
 "within", &
  2*spacing(real(kappa_epsilon))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine kappa_epsilon_matches_epsilonss_per_epsilonee




 subroutine reduces_to_lennard_jones

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
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (lj_potential &
        +2*spacing(real(lj_potential)) ) &
        .ge. &
        (potential(ui, uj, rij)) &
            .and. &
     (lj_potential &
      -2*spacing(real(lj_potential)) ) &
      .le. &
       (potential(ui, uj, rij)) )) then
      print *, " *Assert_Real_Equal failed* in test reduces_to_lennard_jones &
              &[gayberne.fun:78]"
      print *, "  ", "potential(ui, uj, rij) (", &
 potential(ui, uj, rij), &
  ") is not", &
 lj_potential,&
 "within", &
  2*spacing(real(lj_potential))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine reduces_to_lennard_jones




 subroutine separation_millionth_of_tiny_single_precision

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
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(huge(r_absolute) > potential(ui, uj, rij))) then
      print *, " *Assert_True failed* in test separation_millionth_of_tiny_single_precision &
              &[gayberne.fun:102]"
      print *, "  ", "huge(r_absolute) > potential(ui, uj, rij) is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(sqrt(dot_product(rij, rij)) > tiny(r_absolute))) then
      print *, " *Assert_True failed* in test separation_millionth_of_tiny_single_precision &
              &[gayberne.fun:103]"
      print *, "  ", "sqrt(dot_product(rij, rij)) > tiny(r_absolute) is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine separation_millionth_of_tiny_single_precision


 subroutine Setup
  noAssertFailed = .true.
 end subroutine Setup


 subroutine Teardown
 end subroutine Teardown


 subroutine test_gayberne( nTests, nAsserts, nAssertsTested, nFailures )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures

  continue

  call Setup
  call potential_matches_zero_at_cross_contact
  call Teardown

  call Setup
  call kappa_sigma_matches_sigmaee_per_sigmass
  call Teardown

  call Setup
  call kappa_epsilon_matches_epsilonss_per_epsilonee
  call Teardown

  call Setup
  call reduces_to_lennard_jones
  call Teardown

  call Setup
  call separation_millionth_of_tiny_single_precision
  call Teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_gayberne

end module gayberne_fun
