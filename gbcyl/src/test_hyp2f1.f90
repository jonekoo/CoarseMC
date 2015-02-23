program test_hyp2f1
use, intrinsic :: iso_c_binding
use m_hyp2f1_negint, only: hyp2f1_negint
implicit none

interface 
   real(c_double) function hyp2f1(a, b, c, z) bind(C, name='hyp2f1')
     use, intrinsic :: iso_c_binding
     real(c_double), value :: a, b, c, z
   end function hyp2f1
end interface

integer, parameter :: n_m = 14
real(c_double) :: a, b, c, z, y
integer :: n = 100, i, j
integer :: ms(n_m)
ms(1:2) = [4, 10]
ms(3:n_m) = [(j, j = 11, 11 + n_m - 3)]

c = 1
write(*, '("# hyp2f1(-(m - 1) / 2, -(m - 3) / 2, 1, z^2)")')
write(*, '("#",2(' // fmt_char_dp() //', 1X))') '', " m =" // repeat(' ', fmt_char_dp_len() - 2)
write(*, '(' // fmt_char_dp() // ')', advance='no') "# z" // repeat(' ', fmt_char_dp_len() - 2)

do i = 1, n_m
   write(*, '(", "' // fmt_char_dp() // ')', advance='no') ms(i) 
end do
write(*, '(/)', advance='no')

do i = 0, n - 1
  z = i * 1.0 / n
  write(*, '(' // fmt_char_dp() // ')', advance='no') z
  do j = 1, n_m
     a = -(ms(j) - 1) / 2.0
     b = a + 1
     if (mod(ms(j), 2) == 0) then
        y = hyp2f1(a, b, c, z**2)
     else
        !y = hyp2f1_negint(-(ms(j) - 1) / 2, 1, z**2)
        y = hyp2f1_negint(-(ms(j) - 1) / 2, -(ms(j) - 3) / 2, 1, z**2)
     end if
     write(*, '(", "' // fmt_char_dp() // ')', advance='no') y 
  end do
  write(*, '(/)', advance='no')
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
end function

pure function fmt_char_dp_len() result(w)
  integer :: e
  integer :: w
  integer :: d
  real(c_double) :: u 
  e = int(log10(real(range(u)))) + 1
  d = precision(u)
  w = 3 + d + 2 + e
end function

end program
