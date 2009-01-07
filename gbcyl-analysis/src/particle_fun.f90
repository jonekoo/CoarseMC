! particle_fun.f90 - a unit test suite for particle.f90
!
! funit generated this file from particle.fun
! at Wed Dec 10 15:11:23 +0200 2008

module particle_fun

 use particle

 implicit none

 logical :: noAssertFailed

 public :: test_particle

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0



 contains

 subroutine unit_vector_in_origo

  use nrtype, only: dp
  type(particledat) :: particle
  real(dp), dimension(3) :: unit_vector
  particle%x = 0.0 
  particle%y = 0.0
  particle%z = 0.0
  particle%ux = 0.0
  particle%uy = 0.0
  particle%uz = 1.0
  particle%rod = .true.
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%ux &
        +2*spacing(real(particle%ux)) ) &
        .ge. &
        (unit_vector(1)) &
            .and. &
     (particle%ux &
      -2*spacing(real(particle%ux)) ) &
      .le. &
       (unit_vector(1)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_in_origo &
              &[particle.fun:15]"
      print *, "  ", "unit_vector(1) (", &
 unit_vector(1), &
  ") is not", &
 particle%ux,&
 "within", &
  2*spacing(real(particle%ux))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%uy &
        +2*spacing(real(particle%uy)) ) &
        .ge. &
        (unit_vector(2)) &
            .and. &
     (particle%uy &
      -2*spacing(real(particle%uy)) ) &
      .le. &
       (unit_vector(2)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_in_origo &
              &[particle.fun:16]"
      print *, "  ", "unit_vector(2) (", &
 unit_vector(2), &
  ") is not", &
 particle%uy,&
 "within", &
  2*spacing(real(particle%uy))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%uz &
        +2*spacing(real(particle%uz)) ) &
        .ge. &
        (unit_vector(3)) &
            .and. &
     (particle%uz &
      -2*spacing(real(particle%uz)) ) &
      .le. &
       (unit_vector(3)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_in_origo &
              &[particle.fun:17]"
      print *, "  ", "unit_vector(3) (", &
 unit_vector(3), &
  ") is not", &
 particle%uz,&
 "within", &
  2*spacing(real(particle%uz))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  particle%ux = 1.0
  particle%uy = 0.0
  particle%uz = 0.0
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%ux &
        +2*spacing(real(particle%ux)) ) &
        .ge. &
        (unit_vector(1)) &
            .and. &
     (particle%ux &
      -2*spacing(real(particle%ux)) ) &
      .le. &
       (unit_vector(1)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_in_origo &
              &[particle.fun:22]"
      print *, "  ", "unit_vector(1) (", &
 unit_vector(1), &
  ") is not", &
 particle%ux,&
 "within", &
  2*spacing(real(particle%ux))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%uy &
        +2*spacing(real(particle%uy)) ) &
        .ge. &
        (unit_vector(2)) &
            .and. &
     (particle%uy &
      -2*spacing(real(particle%uy)) ) &
      .le. &
       (unit_vector(2)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_in_origo &
              &[particle.fun:23]"
      print *, "  ", "unit_vector(2) (", &
 unit_vector(2), &
  ") is not", &
 particle%uy,&
 "within", &
  2*spacing(real(particle%uy))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%uz &
        +2*spacing(real(particle%uz)) ) &
        .ge. &
        (unit_vector(3)) &
            .and. &
     (particle%uz &
      -2*spacing(real(particle%uz)) ) &
      .le. &
       (unit_vector(3)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_in_origo &
              &[particle.fun:24]"
      print *, "  ", "unit_vector(3) (", &
 unit_vector(3), &
  ") is not", &
 particle%uz,&
 "within", &
  2*spacing(real(particle%uz))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  particle%ux = 1.0/sqrt(2.0_dp)
  particle%uy = 0.0
  particle%uz = 1.0/sqrt(2.0_dp)
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%ux &
        +2*spacing(real(particle%ux)) ) &
        .ge. &
        (unit_vector(1)) &
            .and. &
     (particle%ux &
      -2*spacing(real(particle%ux)) ) &
      .le. &
       (unit_vector(1)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_in_origo &
              &[particle.fun:29]"
      print *, "  ", "unit_vector(1) (", &
 unit_vector(1), &
  ") is not", &
 particle%ux,&
 "within", &
  2*spacing(real(particle%ux))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%uy &
        +2*spacing(real(particle%uy)) ) &
        .ge. &
        (unit_vector(2)) &
            .and. &
     (particle%uy &
      -2*spacing(real(particle%uy)) ) &
      .le. &
       (unit_vector(2)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_in_origo &
              &[particle.fun:30]"
      print *, "  ", "unit_vector(2) (", &
 unit_vector(2), &
  ") is not", &
 particle%uy,&
 "within", &
  2*spacing(real(particle%uy))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%uz &
        +2*spacing(real(particle%uz)) ) &
        .ge. &
        (unit_vector(3)) &
            .and. &
     (particle%uz &
      -2*spacing(real(particle%uz)) ) &
      .le. &
       (unit_vector(3)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_in_origo &
              &[particle.fun:31]"
      print *, "  ", "unit_vector(3) (", &
 unit_vector(3), &
  ") is not", &
 particle%uz,&
 "within", &
  2*spacing(real(particle%uz))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine unit_vector_in_origo




 subroutine unit_vector_x_translation

  type(particledat) :: particle
  real(dp), dimension(3) :: unit_vector
  particle%ux = 0.0
  particle%uy = 1.0
  particle%uz = 0.0
  particle%x = 4.7326
  particle%y = 0.0
  particle%z = 0.0
  particle%rod = .true.
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%ux &
        +2*spacing(real(particle%ux)) ) &
        .ge. &
        (unit_vector(1)) &
            .and. &
     (particle%ux &
      -2*spacing(real(particle%ux)) ) &
      .le. &
       (unit_vector(1)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_x_translation &
              &[particle.fun:47]"
      print *, "  ", "unit_vector(1) (", &
 unit_vector(1), &
  ") is not", &
 particle%ux,&
 "within", &
  2*spacing(real(particle%ux))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%uy &
        +2*spacing(real(particle%uy)) ) &
        .ge. &
        (unit_vector(2)) &
            .and. &
     (particle%uy &
      -2*spacing(real(particle%uy)) ) &
      .le. &
       (unit_vector(2)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_x_translation &
              &[particle.fun:48]"
      print *, "  ", "unit_vector(2) (", &
 unit_vector(2), &
  ") is not", &
 particle%uy,&
 "within", &
  2*spacing(real(particle%uy))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particle%uz &
        +2*spacing(real(particle%uz)) ) &
        .ge. &
        (unit_vector(3)) &
            .and. &
     (particle%uz &
      -2*spacing(real(particle%uz)) ) &
      .le. &
       (unit_vector(3)) )) then
      print *, " *Assert_Real_Equal failed* in test unit_vector_x_translation &
              &[particle.fun:49]"
      print *, "  ", "unit_vector(3) (", &
 unit_vector(3), &
  ") is not", &
 particle%uz,&
 "within", &
  2*spacing(real(particle%uz))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine unit_vector_x_translation




 subroutine unit_vector_y_translation

  type(particledat) :: particle
  real(dp), dimension(3) :: unit_vector
  particle%ux = -1.0
  particle%uy = 0.0
  particle%uz = 0.0
  particle%x = 0.0
  particle%y = 345.8972
  particle%z = 0.0
  particle%rod = .true.
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((unit_vector(1) &
     +1e-7) &
     .ge. &
     (particle%uy) &
             .and. &
     (unit_vector(1) &
     -1e-7) &
     .le. &
     (particle%uy) )) then
      print *, " *Assert_Equal_Within failed* in test unit_vector_y_translation &
              &[particle.fun:65]"
      print *, "  ", "particle%uy (",particle%uy,") is not", &
 unit_vector(1),"within",1e-7
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
    if (.not.((unit_vector(2) &
     +1e-7) &
     .ge. &
     (-particle%ux) &
             .and. &
     (unit_vector(2) &
     -1e-7) &
     .le. &
     (-particle%ux) )) then
      print *, " *Assert_Equal_Within failed* in test unit_vector_y_translation &
              &[particle.fun:66]"
      print *, "  ", "-particle%ux (",-particle%ux,") is not", &
 unit_vector(2),"within",1e-7
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
    if (.not.((unit_vector(3) &
     +1e-7) &
     .ge. &
     (particle%uz) &
             .and. &
     (unit_vector(3) &
     -1e-7) &
     .le. &
     (particle%uz) )) then
      print *, " *Assert_Equal_Within failed* in test unit_vector_y_translation &
              &[particle.fun:67]"
      print *, "  ", "particle%uz (",particle%uz,") is not", &
 unit_vector(3),"within",1e-7
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine unit_vector_y_translation


 subroutine Setup
  noAssertFailed = .true.
 end subroutine Setup


 subroutine Teardown
 end subroutine Teardown


 subroutine test_particle( nTests, nAsserts, nAssertsTested, nFailures )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures

  continue

  call Setup
  call unit_vector_in_origo
  call Teardown

  call Setup
  call unit_vector_x_translation
  call Teardown

  call Setup
  call unit_vector_y_translation
  call Teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_particle

end module particle_fun
