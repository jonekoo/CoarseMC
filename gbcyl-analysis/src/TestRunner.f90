
! TestRunner.f90 - runs fUnit test suites
!
! funit generated this file on Wed Dec 10 18:05:39 +0200 2008.

program TestRunner

    use particle_fun
    use utils_fun
  
  implicit none

  integer, dimension(2) :: numTests, numAsserts, numAssertsTested, numFailures

    write(*,*)
  write(*,*) "particle test suite:"
  call test_particle &
    ( numTests(1), numAsserts(1), numAssertsTested(1), numFailures(1) )
  write(*,*) "Passed", numAssertsTested(1), "of", numAsserts(1), &
             "possible asserts comprising", numTests(1)-numFailures(1), &
             "of", numTests(1), "tests."
    write(*,*)
  write(*,*) "utils test suite:"
  call test_utils &
    ( numTests(2), numAsserts(2), numAssertsTested(2), numFailures(2) )
  write(*,*) "Passed", numAssertsTested(2), "of", numAsserts(2), &
             "possible asserts comprising", numTests(2)-numFailures(2), &
             "of", numTests(2), "tests."
  
  write(*,*)
  write(*,'(a)') "==========[ SUMMARY ]=========="
      write(*,'(a7)',advance="no") " particle:"
  if ( numFailures(1) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a7)',advance="no") " utils:"
  if ( numFailures(2) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,*)

  if ( sum(numFailures) /= 0 ) stop 1

end program TestRunner
