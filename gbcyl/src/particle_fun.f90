! particle_fun.f90 - a unit test suite for particle.f90
!
! funit generated this file from particle.fun
! at Mon Oct 19 15:46:13 +0300 2009

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

 subroutine create

use nrtype
character(len = *), parameter :: particle_string = &
  "gb -0.123e-4 5.677e8 3.232 -0.23423 1e-2 -3e-01"
type(particledat) :: a_particle
real(dp), dimension(3) :: pos
real(dp), dimension(3) :: ori
a_particle = create_particle(particle_string)
pos = position(a_particle)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (-0.123e-4_dp &
        +2*spacing(real(-0.123e-4_dp)) ) &
        .ge. &
        (pos(1)) &
            .and. &
     (-0.123e-4_dp &
      -2*spacing(real(-0.123e-4_dp)) ) &
      .le. &
       (pos(1)) )) then
      print *, " *Assert_Real_Equal failed* in test create &
              &[particle.fun:12]"
      print *, "  ", "pos(1) (", &
 pos(1), &
  ") is not", &
 -0.123e-4_dp,&
 "within", &
  2*spacing(real(-0.123e-4_dp))
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
    if (.not.( (5.677e8_dp &
        +2*spacing(real(5.677e8_dp)) ) &
        .ge. &
        (pos(2)) &
            .and. &
     (5.677e8_dp &
      -2*spacing(real(5.677e8_dp)) ) &
      .le. &
       (pos(2)) )) then
      print *, " *Assert_Real_Equal failed* in test create &
              &[particle.fun:13]"
      print *, "  ", "pos(2) (", &
 pos(2), &
  ") is not", &
 5.677e8_dp,&
 "within", &
  2*spacing(real(5.677e8_dp))
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
    if (.not.( (3.232_dp &
        +2*spacing(real(3.232_dp)) ) &
        .ge. &
        (pos(3)) &
            .and. &
     (3.232_dp &
      -2*spacing(real(3.232_dp)) ) &
      .le. &
       (pos(3)) )) then
      print *, " *Assert_Real_Equal failed* in test create &
              &[particle.fun:14]"
      print *, "  ", "pos(3) (", &
 pos(3), &
  ") is not", &
 3.232_dp,&
 "within", &
  2*spacing(real(3.232_dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
ori = orientation(a_particle)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (-0.23423_dp &
        +2*spacing(real(-0.23423_dp)) ) &
        .ge. &
        (ori(1)) &
            .and. &
     (-0.23423_dp &
      -2*spacing(real(-0.23423_dp)) ) &
      .le. &
       (ori(1)) )) then
      print *, " *Assert_Real_Equal failed* in test create &
              &[particle.fun:16]"
      print *, "  ", "ori(1) (", &
 ori(1), &
  ") is not", &
 -0.23423_dp,&
 "within", &
  2*spacing(real(-0.23423_dp))
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
    if (.not.( (1e-2_dp &
        +2*spacing(real(1e-2_dp)) ) &
        .ge. &
        (ori(2)) &
            .and. &
     (1e-2_dp &
      -2*spacing(real(1e-2_dp)) ) &
      .le. &
       (ori(2)) )) then
      print *, " *Assert_Real_Equal failed* in test create &
              &[particle.fun:17]"
      print *, "  ", "ori(2) (", &
 ori(2), &
  ") is not", &
 1e-2_dp,&
 "within", &
  2*spacing(real(1e-2_dp))
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
    if (.not.( (-3e-01_dp &
        +2*spacing(real(-3e-01_dp)) ) &
        .ge. &
        (ori(3)) &
            .and. &
     (-3e-01_dp &
      -2*spacing(real(-3e-01_dp)) ) &
      .le. &
       (ori(3)) )) then
      print *, " *Assert_Real_Equal failed* in test create &
              &[particle.fun:18]"
      print *, "  ", "ori(3) (", &
 ori(3), &
  ") is not", &
 -3e-01_dp,&
 "within", &
  2*spacing(real(-3e-01_dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine create


 subroutine write

use nrtype
type(particledat) :: particle_written
type(particledat) :: particle_read
character(len = 150) :: particle_string
integer, parameter :: particle_unit = 14
real(dp), dimension(3) :: written_position
real(dp), dimension(3) :: read_position
open(file = 'test_write.tmp', unit = particle_unit, action = 'READWRITE')
particle_written = new_particle()
call write_particle(particle_unit, particle_written) 
rewind(particle_unit)
read(particle_unit, '(A150)') particle_string
particle_read = create_particle(particle_string)
written_position = position(particle_written)
read_position = position(particle_read)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (written_position(1) &
        +2*spacing(real(written_position(1))) ) &
        .ge. &
        (read_position(1)) &
            .and. &
     (written_position(1) &
      -2*spacing(real(written_position(1))) ) &
      .le. &
       (read_position(1)) )) then
      print *, " *Assert_Real_Equal failed* in test write &
              &[particle.fun:37]"
      print *, "  ", "read_position(1) (", &
 read_position(1), &
  ") is not", &
 written_position(1),&
 "within", &
  2*spacing(real(written_position(1)))
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
    if (.not.( (written_position(2) &
        +2*spacing(real(written_position(2))) ) &
        .ge. &
        (read_position(2)) &
            .and. &
     (written_position(2) &
      -2*spacing(real(written_position(2))) ) &
      .le. &
       (read_position(2)) )) then
      print *, " *Assert_Real_Equal failed* in test write &
              &[particle.fun:38]"
      print *, "  ", "read_position(2) (", &
 read_position(2), &
  ") is not", &
 written_position(2),&
 "within", &
  2*spacing(real(written_position(2)))
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
    if (.not.( (written_position(3) &
        +2*spacing(real(written_position(3))) ) &
        .ge. &
        (read_position(3)) &
            .and. &
     (written_position(3) &
      -2*spacing(real(written_position(3))) ) &
      .le. &
       (read_position(3)) )) then
      print *, " *Assert_Real_Equal failed* in test write &
              &[particle.fun:39]"
      print *, "  ", "read_position(3) (", &
 read_position(3), &
  ") is not", &
 written_position(3),&
 "within", &
  2*spacing(real(written_position(3)))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine write


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
  call create
  call Teardown

  call Setup
  call write
  call Teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_particle

end module particle_fun
