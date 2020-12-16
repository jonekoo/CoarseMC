test_suite box

test stringtobox
use utils
use nrtype
character(len = 200) :: boxstring
real(dp) :: lx = 11.1_dp
real(dp) :: ly = 22.2_dp
real(dp) :: lz = 33.3_dp
type(boxdat) :: bp
write(boxstring, '(A, 3' // fmt_char_dp() //', 3L2)') 'box ', lx, ly, lz, .true., .false., .true. 
call createbox(bp, boxstring)
assert_real_equal(getx(bp), lx)
assert_real_equal(gety(bp), ly)
assert_real_equal(getz(bp), lz)
!! :TODO: Test periodicities
end test

end test_suite