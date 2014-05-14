module odd_n_integral
use nrtype
implicit none

contains

function ljwall_tensor(k, radiusA, densityA, epsilonratio, &
  sigmaratio, powers, isotropic_coeffs, anisotropy_coeffs) result(local_tensor)
  real(dp), intent(in) :: k, radiusA, densityA, epsilonratio, sigmaratio
  integer, intent(in) :: powers(:)
  real(dp), intent(in), optional :: isotropic_coeffs(:)
  real(dp), intent(in), optional :: anisotropy_coeffs(:)
  real(dp) :: local_tensor(3, 3) 
  real(dp) :: R
  real(dp) :: isotropic, anisotropy
  integer :: i, n
  R = radiusA / sigmaratio
  local_tensor = 0._dp
  isotropic = 0._dp
  anisotropy = 0._dp
  if (.not. present(isotropic_coeffs) .and. .not. present(anisotropy_coeffs)) &
       return
  do i = 1, size(powers)
    n = powers(i)
    if (present(isotropic_coeffs)) isotropic = isotropic_coeffs(i)
    if (present(anisotropy_coeffs)) anisotropy = anisotropy_coeffs(i)
    
    local_tensor(1, 1) = local_tensor(1, 1) + prefactor(n, R, k) * &
         ((isotropic - anisotropy / 3._dp) * angle_integral(n - 3, 0, k) &
         + real(n - 1, dp) / real(n, dp) * anisotropy * &
         angle_integral(n - 3, 2, k))
    local_tensor(2, 2) = local_tensor(2, 2) + prefactor(n, R, k) * &
         ((isotropic - anisotropy * (1._dp / 3._dp - &
         real(n - 1, dp) / n)) * angle_integral(n - 3, 0, k) - &
         real(n - 1, dp) / n * anisotropy * angle_integral(n - 3, 2, k))  
    local_tensor(3, 3) = local_tensor(3, 3) + prefactor(n, R, k) * &
         ((isotropic - (1._dp / 3._dp - 1._dp / n) * anisotropy) * &
         angle_integral(n - 3, 0, k))
  end do
  local_tensor = epsilonratio * densityA * local_tensor
end function

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
