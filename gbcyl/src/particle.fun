test_suite particle

test create
use nrtype
character(len = *), parameter :: particle_string = &
  "gb -0.123e-4 5.677e8 3.232 -0.23423 1e-2 -3e-01"
type(particledat) :: a_particle
real(dp), dimension(3) :: pos
real(dp), dimension(3) :: ori
a_particle = create_particle(particle_string)
pos = position(a_particle)
assert_real_equal(-0.123e-4_dp, pos(1))
assert_real_equal(5.677e8_dp, pos(2))
assert_real_equal(3.232_dp, pos(3))
ori = orientation(a_particle)
assert_real_equal(-0.23423_dp, ori(1))
assert_real_equal(1e-2_dp, ori(2))
assert_real_equal(-3e-01_dp, ori(3)) 
end test

test write
use nrtype
type(particledat) :: particle_written
type(particledat) :: particle_read
character(len = 150) :: particle_string
integer, parameter :: particle_unit = 14
real(dp), dimension(3) :: written_position
real(dp), dimension(3) :: read_position
open(file = 'test_write.tmp', unit = particle_unit, action = 'READWRITE')
particle_written = new_particle()
call write_particle(particle_unit, particle_written) 
rewind(particle_unit)
read(particle_unit, '(A150)') particle_string
particle_read = create_particle(particle_string)
written_position = position(particle_written)
read_position = position(particle_read)
assert_real_equal(written_position(1), read_position(1))
assert_real_equal(written_position(2), read_position(2))
assert_real_equal(written_position(3), read_position(3))
end test

end test_suite