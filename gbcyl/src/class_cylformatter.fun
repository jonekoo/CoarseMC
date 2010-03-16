test_suite class_cylformatter

test reading
  use particle
  use nrtype
  implicit none
  type(cylformatter) :: cf
  character(len = *), parameter :: statefile = 'testconfiguration.txt'
  type(particledat), dimension(:), pointer :: particles
  real(dp) :: r
  real(dp) :: h
  real(dp), dimension(3) :: pos
  real(dp), dimension(3) :: ori
  integer :: n_particles
  cf = new_cylformatter(statefile)
  call readstate(cf, particles, n_particles, r, h)
  assert_real_equal(9.0_dp, r)
  assert_real_equal(21._dp, h)
  assert_equal(2, n_particles)
  pos = position(particles(1))
  assert_equal(1._dp, pos(1))
  assert_equal(3._dp, pos(2))
  assert_equal(5._dp, pos(3))
  pos = position(particles(2))
  assert_equal(2._dp, pos(1))
  assert_equal(4._dp, pos(2))
  assert_equal(6._dp, pos(3))
  ori = orientation(particles(1))
  assert_equal(0._dp, ori(1))
  assert_equal(0._dp, ori(2))
  assert_equal(1._dp, ori(3))
  ori = orientation(particles(2))
  assert_equal(0._dp, ori(1))
  assert_equal(1._dp, ori(2))
  assert_equal(0._dp, ori(3))  
  call delete(cf)
  !! Test a corrupt configuration
end test

test readlast
  use particle
  use nrtype
  implicit none
  type(cylformatter) :: cf
  character(len = 80) :: statefile = 'testconfiguration.txt'
  type(particledat), dimension(:), pointer :: particles
  real(dp) :: r
  real(dp) :: h
  real(dp), dimension(3) :: ori
  integer :: n_particles
  logical :: isfound
  !! Test with a file that has two configurations of two molecules.
  cf = new_cylformatter(statefile)
  call findlast(cf, isfound)
  assert_true(isfound)
  call readstate(cf, particles, n_particles, r, h)
  ori = orientation(particles(1))
  assert_real_equal(1._dp, ori(2))
  assert_real_equal(0._dp, ori(3))
  ori = orientation(particles(2))
  assert_real_equal(1._dp, ori(3))
  assert_real_equal(0._dp, ori(2))
  !! Test with a file that has a configuration with no beginning
  statefile = 'cylformat-nobegin.txt'
  cf = new_cylformatter(statefile)
  call findlast(cf, isfound)
  assert_false(isfound)
  !! End is there always since it's end-of-file in case of cylformatter.
end test

test generalfindlast
  integer :: readunit = 13
  character(len = 5) :: line
  logical :: isfound
  open(file = 'findlasttest.txt', unit = readunit, action = 'read')
  call findlast(readunit, 'geometry', 'end geometry', isfound)
  assert_true(isfound)
  read(readunit, fmt = '(A5)') line
  assert_equal('geome', line)
  read(readunit, fmt = '(A5)') line
  assert_equal('Jonne', line)
  close(readunit)
  !! Test EOF ending
  open(file = 'findlasttest.txt', unit = readunit, action = 'read')
  call findlast(readunit, 'geometry', 'EOF', isfound)
  assert_true(isfound)
  read(readunit, fmt = '(A5)') line
  assert_equal('geome', line)
  read(readunit, fmt = '(A5)') line
  assert_equal('Jonni', line)
  close(readunit)
end test

end test_suite