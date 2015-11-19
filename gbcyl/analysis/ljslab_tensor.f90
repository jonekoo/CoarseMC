module ljslab_tensor
implicit none

contains

pure function isotropic_integral(n, k, R)
  use num_kind, only: dp
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: k, R
  real(dp) :: isotropic_integral
  real(dp), parameter :: pi = 4. * atan(1.)
  real(dp) :: rw 
  rw = (1. - k) * R
  isotropic_integral = 2. * pi * rw**(3 - n) / (6. - 5 * n + n**2)
end function

pure function anisotropic_integral(n, k, R)
  use num_kind, only: dp
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: k, R
  real(dp) :: anisotropic_integral(3)
  real(dp), parameter :: pi = 4. * atan(1.)
  real(dp) :: rw 
  rw = (1. - k) * R
  anisotropic_integral = [2., -1., -1.] * pi * rw**(3 - n) / (n * (n - 2))
end function

end module ljslab_tensor


pure function ljwall_tensor(k, radiusA, densityA, epsilonratio, sigmaratio, &
     powers, isotropic_coeffs, anisotropy_coeffs) result(local_tensor)
  use ljslab_tensor
  use num_kind, only: dp
  implicit none
  real(dp), intent(in) :: k, radiusA, densityA, epsilonratio, sigmaratio
  integer, intent(in) :: powers(:)
  real(dp), intent(in), optional :: isotropic_coeffs(:)
  real(dp), intent(in), optional :: anisotropy_coeffs(:)
  real(dp) :: local_tensor(3, 3) 
  real(dp) :: R, isotropic, anisotropic(3)
  integer :: i, n
  R = radiusA
  isotropic = 0.
  anisotropic = 0.
  if (present(isotropic_coeffs)) then
     do i = 1, size(powers)
        n = powers(i)
        isotropic = isotropic + sigmaratio**n * isotropic_integral(n, k, R) &
             * isotropic_coeffs(i)
     end do
     do i = 1, 3
        local_tensor(i, i) = isotropic
     end do
  end if
  if (present(anisotropy_coeffs)) then
     do i = 1, size(powers)
        n = powers(i)
        !! Anisotropic integral is the integral of (the anisotropy
        !! times the orientation tensor) over the wall volume.
        anisotropic = anisotropic + sigmaratio**n * &
             anisotropic_integral(n, k, R) * anisotropy_coeffs(i)
     end do
     do i = 1, 3
        local_tensor(i, i) = local_tensor(i, i) + 2. / 3. * anisotropic(i)
     end do
  end if
  local_tensor = local_tensor * epsilonratio * densityA
end function ljwall_tensor

