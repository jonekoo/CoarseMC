module class_poly_box_test
use ftnunit
use class_poly_box
use nrtype
implicit none

contains

subroutine test_all
  call test(test_writeread, "Writing a box to a file and reading it back.") 
end subroutine

subroutine test_writeread
character(len = 200) :: boxstring
real(dp) :: lx = 11.1_dp
real(dp) :: ly = 22.2_dp
real(dp) :: lz = 33.3_dp
real(dp) :: margin = 1.e-9_dp
integer, parameter :: unit = 999
integer :: ios
character(len=80) :: line
type(poly_box) :: boxwritten
type(poly_box) :: boxread
call setx(boxwritten, lx)
call sety(boxwritten, ly)
call setz(boxwritten, lz)
open(unit, status='SCRATCH')
call write(unit, boxwritten)
backspace unit
call read(unit, boxread, ios) 
close(unit)
call assert_equal(0, ios, "Box could not be read")
call assert_comparable(lx, getx(boxread), margin, "Dimensions in x direction &
not comparable")
call assert_comparable(ly, gety(boxread), margin, "Dimensions in y direction &
not comparable")
call assert_comparable(lz, getz(boxread), margin, "Dimensions in z direction &
not comparable")
end subroutine

end module
