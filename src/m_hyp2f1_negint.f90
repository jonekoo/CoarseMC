!> Procedures to compute the Hypergeometric function _2F_1 in the
!! special case in which a or b is a negative integer.
module m_hyp2f1_negint
use iso_c_binding
implicit none

contains

!> Computes the value of the Gauss hypergeometric function in the
!! special case in which @p a or @p b is a negative integer.
!!
!! @see Equation 15.4.1 from Abramowitz and Stegun: Handbook of
!! Mathematical Functions.
!!
pure function hyp2f1_negint(a, b, c, z) result(sum)
  integer, intent(in) :: a, b, c
  real(c_double), intent(in) :: z
  real(c_double) :: sum
  integer(c_int) :: n
  integer :: i

  if (a < 0) n = -a
  if (b < 0 .and. b > a) n = -b

  sum = coeff(a, b, c, n)
  do i = n - 1, 0, -1
     sum = sum * z + coeff(a, b, c, i)
  end do
  
end function

pure function coeff(a, b, c, n) result(res)
  integer, intent(in) :: a, b, c, n
  real(c_double) :: res
  if (n == 0) then
     res = 1.0
  else
     res = ph(a, n) / real(ph(c, n), c_double) * ph(b, n) / &
          real(fact(n), c_double)
  end if
end function coeff

pure recursive function ph(a, n) result(res)
  integer, intent(in) :: a, n
  real(c_double) :: res
  if (n == 0 ) then
     res = 1
  else
     res = a * ph(a + 1, n - 1)
  end if
end function ph

pure recursive function fact(n) result(res)
  integer, intent(in) :: n
  real(c_double) :: res
  if (n > 1) then
     res = n * fact(n - 1)
  else
     res = 1
  end if
end function fact


end module m_hyp2f1_negint
