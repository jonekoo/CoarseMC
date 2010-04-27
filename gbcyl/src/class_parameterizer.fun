test_suite class_parameterizer
!! See file test_parameters.in for test input

test stringparameter
  character(len = *), parameter :: key = 'molfile'
  type(parameterizer) :: p
  character(len = 50) :: value
  p = new_parameterizer('test_parameters.in')
  call getparameter(p, key, value)
  assert_equal(value, 'molecules.in') 
  call delete(p)
end test

test realparameter
  use nrtype
  type(parameterizer) :: p
  character(len = *), parameter :: key = 'gb_kappa'
  real(dp) :: value
  p = new_parameterizer('test_parameters.in')
  call getparameter(p, key, value)
  assert_real_equal(value, 4.4_dp) 
  call delete(p)
end test

test integerparameter
  type(parameterizer) :: p
  character(len = *), parameter :: key = 'mc_sweeps'
  integer :: value
  p = new_parameterizer('test_parameters.in')
  call getparameter(p, key, value)
  assert_equal(value, 1123)
  call delete(p)
end test

end test_suite