module odd_n_integral
use nrtype
implicit none

contains

pure function prefactor(n, R, k) result(res)
  integer, intent(in) :: n
  real(dp), intent(in) :: R, k
  real(dp) :: res
  res = factorial((n - 3) / 2)**2 * 2._dp**(n - 1) / (factorial(n - 2) * &
       (n - 3) * R**(n - 3) * (1 - k**2)**(n - 3))
end function

pure function factorial(n)
  integer, intent(in) :: n
  real(dp) :: factorial
  integer :: i
  factorial = 1._dp
  do i = 1, n
    factorial = factorial * real(i, dp)
  end do
end function

pure function pochhammer(x, n) result(ph)
  real(dp), intent(in) :: x
  integer, intent(in) :: n
  real(dp) :: ph
  integer :: i
  ph = 1
  do i = 0, n - 1
   ph = ph * (x + real(i, dp))
  end do
end function

pure function binomial(n, k) 
  integer, intent(in) :: n, k
  real(dp) :: binomial
  binomial = factorial(n) / (factorial(n - k) * factorial(k))
end function

pure function angle_integral(m, mm, k) result(res)
  integer, intent(in) :: m, mm
  real(dp), intent(in) :: k
  real(dp) :: res
  integer :: l
  real(dp), parameter :: pi = 4._dp * atan(1._dp)
  res = 0._dp
  do l = 0, m, 2
    res = res + binomial(m, l) * k**l * 2._dp**(1 - l - mm) * &
          factorial(l + mm) / factorial((l + mm) / 2)**2 * &
          innersum(m, mm, l, k)
  end do
  res = res * pi / 2._dp
end function

!! This function is probably unnecessary:
pure function angle_integral2(m, mm, k) result(res)
  integer, intent(in) :: m, mm
  real(dp), intent(in) :: k
  real(dp) :: res
  integer :: l
  real(dp), parameter :: pi = 4._dp * atan(1._dp)
  res = 0._dp
  do l = 0, m
    res = res + binomial(m, l) * k**l * (1 + (-1)**l)* 2._dp**(-l - mm) * &
          factorial(l + mm) / factorial((l + mm) / 2)**2 * &
          innersum(m, mm, l, k)
  end do
  res = res * pi / 2._dp
end function

pure function innersum(m, mm, l, k)
    integer, intent(in) :: m, mm, l
    real(dp), intent(in) :: k
    real(dp) :: innersum
    integer :: nn
    innersum = 0._dp
    do nn = 0, (m - l) / 2
      innersum = innersum + pochhammer(real((l - m) / 2, dp), nn) * &
                 pochhammer(0.5_dp, nn) / &
                 pochhammer(real((l + mm + 2) / 2, dp), nn) * k ** (2 * nn) &
                 / factorial(nn)
    end do
end function

!! Restrictions: n has to be odd, mm has to be even!
!!
pure function integral(n, mm, k, R)
  integer, intent(in) :: n
  integer, intent(in) :: mm
  real(dp), intent(in) :: k
  real(dp), intent(in) :: R
  real(dp) :: integral
  integral = 2._dp ** (n - 1) * factorial((n - 3) / 2) ** 2 / &
                   (factorial(n - 2) * (n - 3)) * R ** (3 - n) * &
                   (1 - k ** 2) ** (3 - n) * angle_integral(n - 3, mm, k)
end function


end module
