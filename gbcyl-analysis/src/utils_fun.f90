! utils_fun.f90 - a unit test suite for utils.f90
!
! funit generated this file from utils.fun
! at Wed Dec 10 18:05:39 +0200 2008

module utils_fun

 use utils

 implicit none

 logical :: noAssertFailed

 public :: test_utils

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0





 contains

 subroutine rotation_around_axes

  use nrtype, only: dp
  implicit none
  intrinsic atan
  real(dp), dimension(3) :: ex = (/1.0_dp, 0.0_dp, 0.0_dp/)
  real(dp), dimension(3) :: ey = (/0.0_dp, 1.0_dp, 0.0_dp/)
  real(dp), dimension(3) :: ez = (/0.0_dp, 0.0_dp, 1.0_dp/)
  real(dp), dimension(3) :: u
  real(dp), parameter :: coefficient = -0.988
  real(dp) :: angle
  angle = -2.0*atan(1.0)
  call rotate_vector(coefficient*ex(1), coefficient*ex(2), coefficient*ex(3),&
    & ez(1), ez(2), ez(3), angle, &
    & u(1), u(2), u(3))
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((u(1) &
     +1e-7) &
     .ge. &
     (coefficient*ey(1)) &
             .and. &
     (u(1) &
     -1e-7) &
     .le. &
     (coefficient*ey(1)) )) then
      print *, " *Assert_Equal_Within failed* in test rotation_around_axes &
              &[utils.fun:19]"
      print *, "  ", "coefficient*ey(1) (",coefficient*ey(1),") is not", &
 u(1),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((u(2) &
     +1e-7) &
     .ge. &
     (coefficient*ey(2)) &
             .and. &
     (u(2) &
     -1e-7) &
     .le. &
     (coefficient*ey(2)) )) then
      print *, " *Assert_Equal_Within failed* in test rotation_around_axes &
              &[utils.fun:20]"
      print *, "  ", "coefficient*ey(2) (",coefficient*ey(2),") is not", &
 u(2),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((u(3) &
     +1e-7) &
     .ge. &
     (coefficient*ey(3)) &
             .and. &
     (u(3) &
     -1e-7) &
     .le. &
     (coefficient*ey(3)) )) then
      print *, " *Assert_Equal_Within failed* in test rotation_around_axes &
              &[utils.fun:21]"
      print *, "  ", "coefficient*ey(3) (",coefficient*ey(3),") is not", &
 u(3),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  angle = -angle
  call rotate_vector(ex(1), ex(2), ex(3), ey(1), ey(2), ey(3), angle, &
    & u(1), u(2), u(3))
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((u(1) &
     +1e-7) &
     .ge. &
     (ez(1)) &
             .and. &
     (u(1) &
     -1e-7) &
     .le. &
     (ez(1)) )) then
      print *, " *Assert_Equal_Within failed* in test rotation_around_axes &
              &[utils.fun:25]"
      print *, "  ", "ez(1) (",ez(1),") is not", &
 u(1),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((u(2) &
     +1e-7) &
     .ge. &
     (ez(2)) &
             .and. &
     (u(2) &
     -1e-7) &
     .le. &
     (ez(2)) )) then
      print *, " *Assert_Equal_Within failed* in test rotation_around_axes &
              &[utils.fun:26]"
      print *, "  ", "ez(2) (",ez(2),") is not", &
 u(2),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((u(3) &
     +1e-7) &
     .ge. &
     (ez(3)) &
             .and. &
     (u(3) &
     -1e-7) &
     .le. &
     (ez(3)) )) then
      print *, " *Assert_Equal_Within failed* in test rotation_around_axes &
              &[utils.fun:27]"
      print *, "  ", "ez(3) (",ez(3),") is not", &
 u(3),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  call rotate_vector(ey(1), ey(2), ey(3), ez(1), ez(2), ez(3), angle, &
    & u(1), u(2), u(3))
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((u(1) &
     +1e-7) &
     .ge. &
     (ex(1)) &
             .and. &
     (u(1) &
     -1e-7) &
     .le. &
     (ex(1)) )) then
      print *, " *Assert_Equal_Within failed* in test rotation_around_axes &
              &[utils.fun:30]"
      print *, "  ", "ex(1) (",ex(1),") is not", &
 u(1),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((u(2) &
     +1e-7) &
     .ge. &
     (ex(2)) &
             .and. &
     (u(2) &
     -1e-7) &
     .le. &
     (ex(2)) )) then
      print *, " *Assert_Equal_Within failed* in test rotation_around_axes &
              &[utils.fun:31]"
      print *, "  ", "ex(2) (",ex(2),") is not", &
 u(2),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((u(3) &
     +1e-7) &
     .ge. &
     (ex(3)) &
             .and. &
     (u(3) &
     -1e-7) &
     .le. &
     (ex(3)) )) then
      print *, " *Assert_Equal_Within failed* in test rotation_around_axes &
              &[utils.fun:32]"
      print *, "  ", "ex(3) (",ex(3),") is not", &
 u(3),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine rotation_around_axes


 subroutine Setup
  noAssertFailed = .true.
 end subroutine Setup


 subroutine Teardown
 end subroutine Teardown


 subroutine test_utils( nTests, nAsserts, nAssertsTested, nFailures )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures

  continue

  call Setup
  call rotation_around_axes
  call Teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_utils

end module utils_fun
