program main
   use pFUnit
   use mpi
   use gayberne_pfunit
   use class_simplelist_pfunit
   use cell_energy_pfunit
   implicit none

   external testZeroVectors
   external testSameVector
   external testOrthogonalVectors
   external testFailOnPurpose

  
   type (TestSuite_type) :: suite
   type (TestSuite_type) :: gayberne_suite
   type (TestSuite_type) :: class_simplelist_suite
   type (TestSuite_type) :: cell_energy_suite
   type (TestResult_type) :: result
   character(len=100) :: summary_statement

   call pFUnit_init()

! Build suite from test procedures:
   suite = TestSuite('vector tests')
   gayberne_suite = TestSuite('Gay-Berne potential tests')
   class_simplelist_suite = TestSuite('Simple cell list tests')
   cell_energy_suite = TestSuite('Cell list tests')

   call add(suite, TestCase1Step('testZeroVectors', testZeroVectors))
   call add(suite, TestCase1Step('testSameVector', testSameVector))
   call add(suite, TestCase1Step('testOrthogonalVectors', testOrthogonalVectors))
 !  call add(suite, TestCase1Step('testFailOnPurpose', testFailOnPurpose))
   call add(gayberne_suite, TestCase1Step('Calculate potential at cross contact', test_zero_at_cross_contact))
   call add(class_simplelist_suite, TestCase1Step('Test creation of list',test_new_simplelist))
   call add(class_simplelist_suite, TestCase1Step('Test updating of list',test_updateall))
   call add(class_simplelist_suite, TestCase1Step('Test updating of list one particle at a time', test_updatesingle))
   call add(class_simplelist_suite, TestCase1Step('Test neighbour mask',test_nbrmask))
   call add(cell_energy_suite, TestCase1Step('Test neighbour mask',test_ce_nbrmask))
! Run the tests and accumulate the results in "result"
   result = newTestResult(mode=MODE_USE_STDOUT)
   call Run(suite, result)
   summary_statement=Summary(result)
   print*,trim(summary_statement)

!   call clean(result)
   call clean(suite)

   call Run(gayberne_suite, result)
   summary_statement=Summary(result)
   print*,trim(summary_statement)

!   call clean(result)
   call clean(gayberne_suite)

   call Run(class_simplelist_suite, result)
   summary_statement=Summary(result)
   print*,trim(summary_statement)

!   call clean(result)
   call clean(class_simplelist_suite)

   call Run(cell_energy_suite, result)
   summary_statement=Summary(result)
   print*,trim(summary_statement)

   call clean(result)
   call clean(cell_energy_suite)

   call pFUnit_finalize()

end program main
