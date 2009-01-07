test_suite particle

test unit_vector_in_origo
  use nrtype, only: dp
  type(particledat) :: particle
  real(dp), dimension(3) :: unit_vector
  particle%x = 0.0 
  particle%y = 0.0
  particle%z = 0.0
  particle%ux = 0.0
  particle%uy = 0.0
  particle%uz = 1.0
  particle%rod = .true.
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  assert_real_equal(particle%ux, unit_vector(1))
  assert_real_equal(particle%uy, unit_vector(2))
  assert_real_equal(particle%uz, unit_vector(3))
  particle%ux = 1.0
  particle%uy = 0.0
  particle%uz = 0.0
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  assert_real_equal(particle%ux, unit_vector(1))
  assert_real_equal(particle%uy, unit_vector(2))
  assert_real_equal(particle%uz, unit_vector(3))
  particle%ux = 1.0/sqrt(2.0_dp)
  particle%uy = 0.0
  particle%uz = 1.0/sqrt(2.0_dp)
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  assert_real_equal(particle%ux, unit_vector(1))
  assert_real_equal(particle%uy, unit_vector(2))
  assert_real_equal(particle%uz, unit_vector(3))
end test



test unit_vector_x_translation
  type(particledat) :: particle
  real(dp), dimension(3) :: unit_vector
  particle%ux = 0.0
  particle%uy = 1.0
  particle%uz = 0.0
  particle%x = 4.7326
  particle%y = 0.0
  particle%z = 0.0
  particle%rod = .true.
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  assert_real_equal(particle%ux, unit_vector(1))
  assert_real_equal(particle%uy, unit_vector(2))
  assert_real_equal(particle%uz, unit_vector(3))
end test



test unit_vector_y_translation
  type(particledat) :: particle
  real(dp), dimension(3) :: unit_vector
  particle%ux = -1.0
  particle%uy = 0.0
  particle%uz = 0.0
  particle%x = 0.0
  particle%y = 345.8972
  particle%z = 0.0
  particle%rod = .true.
  call unitvec(particle, unit_vector(1), unit_vector(2), unit_vector(3))
  assert_equal_within(particle%uy, unit_vector(1), 1e-7)
  assert_equal_within(-particle%ux, unit_vector(2), 1e-7)
  assert_equal_within(particle%uz, unit_vector(3), 1e-7)
end test

end test_suite