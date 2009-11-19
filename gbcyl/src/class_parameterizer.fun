test_suite class_parameterizer
!! See file test_parameters.in for test input

test string_parameter
  character(len = *), parameter :: key = 'molfile'
  type(parameterizer) :: p
  character(len = 50) :: value
  p = new_parameterizer('test_parameters.in')
  call get_parameter(p, key, value)
  assert_equal(value, 'molecules.in') 
end test

test real_parameter
  use nrtype
  type(parameterizer) :: p
  character(len = *), parameter :: key = 'gb_kappa'
  real(dp) :: value
  p = new_parameterizer('test_parameters.in')
  call get_parameter(p, key, value)
  assert_real_equal(value, 4.4_dp) 
end test

test integer_parameter
  type(parameterizer) :: p
  character(len = *), parameter :: key = 'mc_sweeps'
  integer :: value
  p = new_parameterizer('test_parameters.in')
  call get_parameter(p, key, value)
  assert_equal(value, 1123)
end test

!test parameter_not_found
!end test

end test_suite