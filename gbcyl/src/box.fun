test_suite box

test string_to_box
use utils
use nrtype
character(len = 200) :: box_string
real(dp) :: lx = 11.1_dp
real(dp) :: ly = 22.2_dp
real(dp) :: lz = 33.3_dp
type(boxdat) :: bp
write(box_string, '(A, 3' // fmt_char_dp() //', 3L2)') 'box ', lx, ly, lz, .true., .false., .true. 
call create_box(bp, box_string)
assert_real_equal(get_x(bp), lx)
assert_real_equal(get_y(bp), ly)
assert_real_equal(get_z(bp), lz)
!! :TODO: Test periodicities
end test

end test_suite