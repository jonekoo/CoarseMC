module utils_test
use ftnunit
use utils
implicit none

contains

subroutine test_all()
  call test(test_substrcount, "Test that substrings are counted correctly.")
  call test(test_splitstr, "Test that string is splitted correctly.")
end subroutine

subroutine test_substrcount()
  character(len=*), parameter :: str = 'osetspojsetpojset'
  call assert_equal(substrcount(str, ''), 0, "More than zero empty substrings found!")
  call assert_equal(substrcount(str, 'o'), 3, "Wrong amount of one char substrings!")
  call assert_equal(substrcount(str, 'set'), 3, "Wrong amount of 3-char substrings!")
end subroutine

subroutine test_splitstr()
  character(len=*), parameter :: str ='xy,xz,y,z,x,'
  character(len=len(str)), dimension(:), allocatable :: splitted
  call splitstr(str, ',', splitted)
  call assert_equal(size(splitted), 6, "Splitted string has wrong length!")
  call assert_true(splitted(1) == 'xy', "Splitted string has wrong first substring!")
  call assert_true(splitted(6) == '', "Splitted string has wrong last substring!")
  deallocate(splitted)
  call splitstr(str, 'z,', splitted)
  call assert_equal(size(splitted), 3, "Splitted string has wrong length!")
  call assert_true(splitted(1) == str(1:4), "Splitted string has wrong first substring!")
  call assert_true(splitted(2) == str(7:8), "Splitted string has wrong second substring!")
  call assert_true(splitted(3) == str(11:12), "Splitted string has wrong third substring!")
end subroutine
 
subroutine test_joinstr()
  character(len=*), parameter :: str_array(5) =(/"xy","xz","y ","z ","x "/)
  character(len=200) :: joined
  call join(str_array, ',', joined)
  call assert_true(trim(adjustl(joined)) == 'xy,xz,y,z,x', "Joined string does not match.")
end subroutine
 
end module


