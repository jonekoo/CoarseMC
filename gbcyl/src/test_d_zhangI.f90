program test_d_zhangI
use, intrinsic :: iso_c_binding
use cylinder_integrals, only: zhangI
implicit none

integer, parameter :: m_lb = 11
integer, parameter :: m_ub = 35
real(c_double) :: y
integer :: j
real(c_double) :: k
integer :: m
real(c_double) :: R = 9.0

do m = m_lb, m_ub
  do j = 0, 19
     k = 0.9 + j * 0.005
     y = d_zhangI(m, k, R)
     write(*, '(4(' // fmt_char_dp() // '))') m, k, R, y 
  end do
end do

contains 


!> Returns a formatting character for a double precision real number. To be 
!! used when consistent formatting of real numbers is needed. 
!!
pure function fmt_char_dp() result(format_char)
  integer :: e
  integer :: w
  integer :: d
  real(c_double) :: u 
  character(len = 50) :: w_char !! Width of field
  character(len = 50) :: d_char !! Width of decimal field
  character(len = 50) :: e_char !! width of exponent
  character(len = 50) :: format_char
  e = int(log10(real(range(u)))) + 1
  d = precision(u)
  w = 3 + d + 2 + e
  write(w_char, *) w
  write(e_char, *) e
  write(d_char, *) d
  write(format_char, *) 'G' // trim(adjustl(w_char)) // '.' // trim(adjustl(d_char)) // 'E' // trim(adjustl(e_char)) 
  format_char = trim(adjustl(format_char))
end function fmt_char_dp

pure function fmt_char_dp_len() result(w)
  integer :: e
  integer :: w
  integer :: d
  real(c_double) :: u 
  e = int(log10(real(range(u)))) + 1
  d = precision(u)
  w = 3 + d + 2 + e
end function fmt_char_dp_len

end program test_zhangI
