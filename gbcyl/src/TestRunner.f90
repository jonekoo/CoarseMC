
! TestRunner.f90 - runs fUnit test suites
!
! funit generated this file on Wed Feb 11 16:42:55 +0200 2009.

program TestRunner

    use gayberne_fun
  
  implicit none

  integer, dimension(1) :: numTests, numAsserts, numAssertsTested, numFailures

    write(*,*)
  write(*,*) "gayberne test suite:"
  call test_gayberne &
    ( numTests(1), numAsserts(1), numAssertsTested(1), numFailures(1) )
  write(*,*) "Passed", numAssertsTested(1), "of", numAsserts(1), &
             "possible asserts comprising", numTests(1)-numFailures(1), &
             "of", numTests(1), "tests."
  
  write(*,*)
  write(*,'(a)') "==========[ SUMMARY ]=========="
      write(*,'(a10)',advance="no") " gayberne:"
  if ( numFailures(1) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,*)

  if ( sum(numFailures) /= 0 ) stop 1

end program TestRunner
