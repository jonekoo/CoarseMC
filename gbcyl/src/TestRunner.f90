
! TestRunner.f90 - runs fUnit test suites
!
! funit generated this file on Tue Mar 16 12:24:02 +0200 2010.

program TestRunner

    use box_fun
    use cell_fun
    use cell_energy_fun
    use class_cylformatter_fun
    use class_parameterizer_fun
    use cylinder_fun
    use gayberne_fun
    use io_fun
    use mtmod_fun
    use particle_fun
    use pt_fun
    use verlet_fun
  
  implicit none

  integer, dimension(12) :: numTests, numAsserts, numAssertsTested, numFailures

    write(*,*)
  write(*,*) "box test suite:"
  call test_box &
    ( numTests(1), numAsserts(1), numAssertsTested(1), numFailures(1) )
  write(*,*) "Passed", numAssertsTested(1), "of", numAsserts(1), &
             "possible asserts comprising", numTests(1)-numFailures(1), &
             "of", numTests(1), "tests."
    write(*,*)
  write(*,*) "cell test suite:"
  call test_cell &
    ( numTests(2), numAsserts(2), numAssertsTested(2), numFailures(2) )
  write(*,*) "Passed", numAssertsTested(2), "of", numAsserts(2), &
             "possible asserts comprising", numTests(2)-numFailures(2), &
             "of", numTests(2), "tests."
    write(*,*)
  write(*,*) "cell_energy test suite:"
  call test_cell_energy &
    ( numTests(3), numAsserts(3), numAssertsTested(3), numFailures(3) )
  write(*,*) "Passed", numAssertsTested(3), "of", numAsserts(3), &
             "possible asserts comprising", numTests(3)-numFailures(3), &
             "of", numTests(3), "tests."
    write(*,*)
  write(*,*) "class_cylformatter test suite:"
  call test_class_cylformatter &
    ( numTests(4), numAsserts(4), numAssertsTested(4), numFailures(4) )
  write(*,*) "Passed", numAssertsTested(4), "of", numAsserts(4), &
             "possible asserts comprising", numTests(4)-numFailures(4), &
             "of", numTests(4), "tests."
    write(*,*)
  write(*,*) "class_parameterizer test suite:"
  call test_class_parameterizer &
    ( numTests(5), numAsserts(5), numAssertsTested(5), numFailures(5) )
  write(*,*) "Passed", numAssertsTested(5), "of", numAsserts(5), &
             "possible asserts comprising", numTests(5)-numFailures(5), &
             "of", numTests(5), "tests."
    write(*,*)
  write(*,*) "cylinder test suite:"
  call test_cylinder &
    ( numTests(6), numAsserts(6), numAssertsTested(6), numFailures(6) )
  write(*,*) "Passed", numAssertsTested(6), "of", numAsserts(6), &
             "possible asserts comprising", numTests(6)-numFailures(6), &
             "of", numTests(6), "tests."
    write(*,*)
  write(*,*) "gayberne test suite:"
  call test_gayberne &
    ( numTests(7), numAsserts(7), numAssertsTested(7), numFailures(7) )
  write(*,*) "Passed", numAssertsTested(7), "of", numAsserts(7), &
             "possible asserts comprising", numTests(7)-numFailures(7), &
             "of", numTests(7), "tests."
    write(*,*)
  write(*,*) "io test suite:"
  call test_io &
    ( numTests(8), numAsserts(8), numAssertsTested(8), numFailures(8) )
  write(*,*) "Passed", numAssertsTested(8), "of", numAsserts(8), &
             "possible asserts comprising", numTests(8)-numFailures(8), &
             "of", numTests(8), "tests."
    write(*,*)
  write(*,*) "mtmod test suite:"
  call test_mtmod &
    ( numTests(9), numAsserts(9), numAssertsTested(9), numFailures(9) )
  write(*,*) "Passed", numAssertsTested(9), "of", numAsserts(9), &
             "possible asserts comprising", numTests(9)-numFailures(9), &
             "of", numTests(9), "tests."
    write(*,*)
  write(*,*) "particle test suite:"
  call test_particle &
    ( numTests(10), numAsserts(10), numAssertsTested(10), numFailures(10) )
  write(*,*) "Passed", numAssertsTested(10), "of", numAsserts(10), &
             "possible asserts comprising", numTests(10)-numFailures(10), &
             "of", numTests(10), "tests."
    write(*,*)
  write(*,*) "pt test suite:"
  call test_pt &
    ( numTests(11), numAsserts(11), numAssertsTested(11), numFailures(11) )
  write(*,*) "Passed", numAssertsTested(11), "of", numAsserts(11), &
             "possible asserts comprising", numTests(11)-numFailures(11), &
             "of", numTests(11), "tests."
    write(*,*)
  write(*,*) "verlet test suite:"
  call test_verlet &
    ( numTests(12), numAsserts(12), numAssertsTested(12), numFailures(12) )
  write(*,*) "Passed", numAssertsTested(12), "of", numAsserts(12), &
             "possible asserts comprising", numTests(12)-numFailures(12), &
             "of", numTests(12), "tests."
  
  write(*,*)
  write(*,'(a)') "==========[ SUMMARY ]=========="
      write(*,'(a8)',advance="no") " box:"
  if ( numFailures(1) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " cell:"
  if ( numFailures(2) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " cell_energy:"
  if ( numFailures(3) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " class_cylformatter:"
  if ( numFailures(4) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " class_parameterizer:"
  if ( numFailures(5) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " cylinder:"
  if ( numFailures(6) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " gayberne:"
  if ( numFailures(7) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " io:"
  if ( numFailures(8) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " mtmod:"
  if ( numFailures(9) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " particle:"
  if ( numFailures(10) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " pt:"
  if ( numFailures(11) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a8)',advance="no") " verlet:"
  if ( numFailures(12) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,*)

  if ( sum(numFailures) /= 0 ) stop 1

end program TestRunner
