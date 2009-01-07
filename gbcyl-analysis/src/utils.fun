test_suite utils



test rotation_around_axes
  use nrtype, only: dp
  implicit none
  intrinsic atan
  real(dp), dimension(3) :: ex = (/1.0_dp, 0.0_dp, 0.0_dp/)
  real(dp), dimension(3) :: ey = (/0.0_dp, 1.0_dp, 0.0_dp/)
  real(dp), dimension(3) :: ez = (/0.0_dp, 0.0_dp, 1.0_dp/)
  real(dp), dimension(3) :: u
  real(dp), parameter :: coefficient = -0.988
  real(dp) :: angle
  angle = -2.0*atan(1.0)
  call rotate_vector(coefficient*ex(1), coefficient*ex(2), coefficient*ex(3),&
    & ez(1), ez(2), ez(3), angle, &
    & u(1), u(2), u(3))
  assert_equal_within(coefficient*ey(1), u(1), 1e-7)
  assert_equal_within(coefficient*ey(2), u(2), 1e-7)
  assert_equal_within(coefficient*ey(3), u(3), 1e-7)
  angle = -angle
  call rotate_vector(ex(1), ex(2), ex(3), ey(1), ey(2), ey(3), angle, &
    & u(1), u(2), u(3))
  assert_equal_within(ez(1), u(1), 1e-7)
  assert_equal_within(ez(2), u(2), 1e-7)
  assert_equal_within(ez(3), u(3), 1e-7)
  call rotate_vector(ey(1), ey(2), ey(3), ez(1), ez(2), ez(3), angle, &
    & u(1), u(2), u(3))
  assert_equal_within(ex(1), u(1), 1e-7)
  assert_equal_within(ex(2), u(2), 1e-7)
  assert_equal_within(ex(3), u(3), 1e-7)
end test




end test_suite