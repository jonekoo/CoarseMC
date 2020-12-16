test_suite particle

test create
use nrtype
character(len = *), parameter :: particlestring = &
  "gb -0.123e-4 5.677e8 3.232 -0.23423 1e-2 -3e-01"
type(particledat) :: aparticle
real(dp), dimension(3) :: pos
real(dp), dimension(3) :: ori
aparticle = createparticle(particlestring)
pos = position(aparticle)
assert_real_equal(-0.123e-4_dp, pos(1))
assert_real_equal(5.677e8_dp, pos(2))
assert_real_equal(3.232_dp, pos(3))
ori = orientation(aparticle)
assert_real_equal(-0.23423_dp, ori(1))
assert_real_equal(1e-2_dp, ori(2))
assert_real_equal(-3e-01_dp, ori(3)) 
end test

test write
use nrtype
type(particledat) :: particlewritten
type(particledat) :: particleread
character(len = 150) :: particlestring
integer, parameter :: particleunit = 14
real(dp), dimension(3) :: writtenposition
real(dp), dimension(3) :: readposition
open(file = 'test_write.tmp', unit = particleunit, action = 'READWRITE')
particlewritten = new_particle()
call writeparticle(particlewritten, particleunit) 
rewind(particleunit)
read(particleunit, '(A150)') particlestring
particleread = createparticle(particlestring)
writtenposition = position(particlewritten)
readposition = position(particleread)
assert_real_equal(writtenposition(1), readposition(1))
assert_real_equal(writtenposition(2), readposition(2))
assert_real_equal(writtenposition(3), readposition(3))
end test

end test_suite