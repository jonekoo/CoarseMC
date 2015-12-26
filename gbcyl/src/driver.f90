program main
   use pFUnit
   !use mpi
   use gayberne_pfunit
   use class_simplelist_pfunit
   use utils_pfunit
   use particle_pfunit
   implicit none

   type (TestSuite_type) :: gayberne_suite
   type (TestSuite_type) :: class_simplelist_suite
   type (TestSuite_type) :: utils_suite
   type (TestSuite_type) :: particle_suite
   type (TestResult_type) :: result
   character(len=100) :: summary_statement

   call pFUnit_init()

   ! Build suite from test procedures:
   gayberne_suite = TestSuite('Gay-Berne potential tests')
   call add(gayberne_suite, TestCase1Step(&
        'Calculate potential at cross contact', test_zero_at_cross_contact))
   call add(gayberne_suite, TestCase1Step(&
        'Check kappa for sigma', test_kappa_sigma_defined))
   call add(gayberne_suite, TestCase1Step(&
        'Check kappa for epsilon', test_kappa_epsilon_defined))
   call add(gayberne_suite, TestCase1Step(&
        'Test GB potential with spherical particles', &
        test_reduces_to_lennard_jones))
   call add(gayberne_suite, TestCase1Step(&
        'Test derivative zeros at different configurations',&
        test_force_zeros))
   call add(gayberne_suite, TestCase1Step('Test small separation',&
        test_small_separation))
   call add(gayberne_suite, TestCase1Step('Test normal operation',&
        test_normalop))
   call add(gayberne_suite, TestCase1Step('Test force vs finite difference',&
        test_forcevsfinitedifference))

   class_simplelist_suite = TestSuite('Simple cell list tests')
   call add(class_simplelist_suite, TestCase1Step('Test creation of list',&
        test_new_simplelist))
   call add(class_simplelist_suite, TestCase1Step('Test updating of list',&
        test_updateall))
   call add(class_simplelist_suite, TestCase1Step('Test neighbour mask',&
        test_nbrmask))

   utils_suite = TestSuite('Tests for utility routines in module utils.')
   call add(utils_suite, TestCase1Step('Test reading a file to a string.', &
        test_readstr))
   !call add(utils_suite, TestCase1Step('Test parsing a file to lines.', &
   !     test_readlines))

   particle_suite = TestSuite('Tests for particle i/o')
   call add(particle_suite, &
        TestCase1Step('Test writing and reading particledat from json.', &
        test_particledat_json_io))
   call add(particle_suite, &
        TestCase1Step('Test writing and reading rod from json.', &
        test_rod_json_io))
   call add(particle_suite, &
        TestCase1Step('Test writing and reading point from json.', &
        test_point_json_io))
   call add(particle_suite, &
        TestCase1Step('Test writing and reading rod array from json.', &
        test_rodarray_json_io))
   call add(particle_suite, &
        TestCase1Step('Test writing and reading point array from json.', &
        test_pointarray_json_io))
   
   
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

   call Run(utils_suite, result)
   summary_statement=Summary(result)
   print*,trim(summary_statement)
   call clean(utils_suite)

   call Run(particle_suite, result)
   summary_statement=Summary(result)
   print*,trim(summary_statement)
   call clean(particle_suite)

   call pFUnit_finalize()

end program main
