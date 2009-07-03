
! TestRunner.f90 - runs fUnit test suites
!
! funit generated this file on Fri Jul 03 16:13:49 +0300 2009.

program TestRunner

    use gayberne_fun
    use pt_fun
  
  implicit none

  integer, dimension(2) :: numTests, numAsserts, numAssertsTested, numFailures

    write(*,*)
  write(*,*) "gayberne test suite:"
  call test_gayberne &
    ( numTests(1), numAsserts(1), numAssertsTested(1), numFailures(1) )
  write(*,*) "Passed", numAssertsTested(1), "of", numAsserts(1), &
             "possible asserts comprising", numTests(1)-numFailures(1), &
             "of", numTests(1), "tests."
    write(*,*)
  write(*,*) "pt test suite:"
  call test_pt &
    ( numTests(2), numAsserts(2), numAssertsTested(2), numFailures(2) )
  write(*,*) "Passed", numAssertsTested(2), "of", numAsserts(2), &
             "possible asserts comprising", numTests(2)-numFailures(2), &
             "of", numTests(2), "tests."
  
  write(*,*)
  write(*,'(a)') "==========[ SUMMARY ]=========="
      write(*,'(a4)',advance="no") " gayberne:"
  if ( numFailures(1) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a4)',advance="no") " pt:"
  if ( numFailures(2) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,*)

  if ( sum(numFailures) /= 0 ) stop 1

end program TestRunner
