function ljwall_tensor(k, radiusA, densityA, epsilonratio, &
  sigmaratio, powers, isotropic_coeffs, anisotropy_coeffs) result(local_tensor)
  use cylinder_integrals, only: zhangI, gamma
  use nrtype, only: dp
  implicit none
  real(dp), intent(in) :: k, radiusA, densityA, epsilonratio, sigmaratio
  integer, intent(in) :: powers(:)
  real(dp), intent(in), optional :: isotropic_coeffs(:)
  real(dp), intent(in), optional :: anisotropy_coeffs(:)
  real(dp) :: local_tensor(3, 3) 
  real(dp) :: R
  real(dp) :: isotropic, anisotropy
  integer :: i, n
  real(dp), parameter :: pi = 4 * atan(1._dp)
  R = radiusA
  local_tensor = 0._dp
  isotropic = 0._dp
  anisotropy = 0._dp
  if (.not. present(isotropic_coeffs) .and. .not. present(anisotropy_coeffs)) &
       return
  do i = 1, size(powers)
    n = powers(i)
    if (present(isotropic_coeffs)) isotropic = isotropic_coeffs(i)
    if (present(anisotropy_coeffs)) anisotropy = anisotropy_coeffs(i)
    
    local_tensor(1, 1) = local_tensor(1, 1) + sigmaratio**n * sqrt(pi) * &
         gamma((n - 1) / 2._dp) / gamma(n / 2._dp) * &
         !gammaxper2(n - 1) / gammaxper2(n) * &
         ((isotropic - anisotropy / 3) * zhangI(n - 2, k, R) + &
         anisotropy * (n - 1) / real(n * (n - 3), dp) * (&
         (k**2 - 1)**2 * R**2 * (n - 1) / (4 * k**2) * zhangI(n, k, R) + &
         (n - 5) / (4 * k**2 * R**2) * zhangI(n - 4, k, R) + &
         (k**2 - 1) * (n - 3) / (2 * k**2) * zhangI(n - 2, k, R)))
    local_tensor(2, 2) = local_tensor(2, 2) + sigmaratio**n * sqrt(pi) * &
         gamma((n - 1) / 2._dp) / gamma(n / 2._dp) * &
         !gammaxper2(n - 1) / gammaxper2(n) * &
         ((isotropic + anisotropy * (2 * n - 3) / (3._dp * n)) * &
         zhangI(n - 2, k, R) - anisotropy * (n - 1) / real(n * (n - 3), dp) &
         * ((k**2 - 1)**2 * R**2 * (n - 1) / (4 * k**2) * zhangI(n, k, R) + &
         (n - 5) / (4 * k**2 * R**2) * zhangI(n - 4, k, R) + &
         (k**2 - 1) * (n - 3) / (2 * k**2) * zhangI(n - 2, k, R)))
    local_tensor(3, 3) = local_tensor(3, 3) + sigmaratio**n * sqrt(pi) * &
         gamma((n - 1) / 2._dp) / gamma(n / 2._dp) * & 
         !gammaxper2(n - 1) / gammaxper2(n) * &
         (isotropic - anisotropy * (n - 3) / (3._dp * n)) * &
         zhangI(n - 2, k, R)
  end do
  local_tensor = epsilonratio * densityA * local_tensor
  contains
    !> Computes gamma(n / 2) using formula 6.1.12 from Abramowitz & 
    !! Stegun, Handbook of Mathematical Functions.
    real(dp) function gammaxper2(n)
      !use cylinder_integrals, only: gamma
      use nrtype, only: dp
      implicit none
      integer, intent(in) :: n
      integer :: j
      gammaxper2 = 1.0
      if (mod(n, 2) == 0) then 
         !! n is even, just compute factorial for n / 2 - 1.
         do j = 2, n / 2 - 1
            gammaxper2 = gammaxper2 * real(j, dp)
         end do
      else
         !! use 6.1.12
         !! gamma(m + 1/2) = gamma(n / 2) = gamma((n-1)/2 + 1/2) => m = (n - 1) / 2 
         !! 2 * m - 1 = n - 2
         do j = 3, n - 2, 2
            gammaxper2 = gammaxper2 * real(j, dp)
         end do
         !! gamma(1/2) = sqrt(pi)
         gammaxper2 = gammaxper2 * sqrt(pi) / 2.0**((n - 1) / 2)
      end if
    end function gammaxper2
end function


