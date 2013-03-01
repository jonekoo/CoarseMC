program main
   use pFUnit
   use mpi
   use gayberne_pfunit
   use class_simplelist_pfunit
   implicit none

   type (TestSuite_type) :: gayberne_suite
   type (TestSuite_type) :: class_simplelist_suite
   type (TestResult_type) :: result
   character(len=100) :: summary_statement

   call pFUnit_init()

! Build suite from test procedures:
   gayberne_suite = TestSuite('Gay-Berne potential tests')
   class_simplelist_suite = TestSuite('Simple cell list tests')

   call add(gayberne_suite, TestCase1Step('Calculate potential at cross contact', test_zero_at_cross_contact))
   call add(gayberne_suite, TestCase1Step('Check kappa for sigma', test_kappa_sigma_defined))
   call add(gayberne_suite, TestCase1Step('Check kappa for epsilon', test_kappa_epsilon_defined))
   call add(gayberne_suite, TestCase1Step('Test GB potential with spherical particles', test_reduces_to_lennard_jones))
   call add(gayberne_suite, TestCase1Step('Test derivative zeros at different configurations', test_derivative_zeros))
   call add(gayberne_suite, TestCase1Step('Test small separation', test_small_separation))
   call add(gayberne_suite, TestCase1Step('Test normal operation', test_normalop))
   call add(gayberne_suite, TestCase1Step('Test force vs finite difference', test_forcevsfinitedifference))
   call add(class_simplelist_suite, TestCase1Step('Test creation of list',test_new_simplelist))
   call add(class_simplelist_suite, TestCase1Step('Test updating of list',test_updateall))
   call add(class_simplelist_suite, TestCase1Step('Test updating of list one particle at a time', test_updatesingle))
   call add(class_simplelist_suite, TestCase1Step('Test neighbour mask',test_nbrmask))

! Run the tests and accumulate the results in "result"
   result = newTestResult(mode=MODE_USE_STDOUT)
   call Run(gayberne_suite, result)
   summary_statement=Summary(result)
   print*,trim(summary_statement)

   call clean(gayberne_suite)

   call Run(class_simplelist_suite, result)
   summary_statement=Summary(result)
   print*,trim(summary_statement)

   call clean(class_simplelist_suite)

   call pFUnit_finalize()

end program main
