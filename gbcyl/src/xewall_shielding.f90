!> Computes the nuclear magnetic shielding tensor for a 129/131Xe atom
!! inside a cylindrical cavity with radius @p radiusA. Cavity wall is
!! infinitely thick and consists of smoothly and evenly distributed
!! Lennard-Jones particles. The local coordinate system has z as the
!! cylinder axis and x pointing to the radial direction. y-axis is
!! perpendicular to both. The @p epsilonratio is used to scale the strength
!! of the shielding with respect to the Xe-Xe shielding. @p sigmaratio is
!! used to scale the distance to the / curvature of the wall compared to a
!! wall consisting of Xe atoms. 
!! 
!! @param k the ratio (distance of Xe from cylinder axis) / 
!!    (inner radius of the cylindrical cavity)
!! @param radiusA the inner radius of the cavity in aengstroms
!! @param densityA the density of the wall in 1/aengstroms^3
!! @param epsilonratio = xewall_epsilon_0 / xexe_epsilon_0
!! @param sigmaratio = xewall_sigma_0 / xexe_sigma_0
!!
!! @return the shielding tensor (ppm) in the local coordinate system. 
!!
function xewall_shielding(k, radiusA, densityA, epsilonratio, &
  sigmaratio) result(local_tensor)
  use num_kind, only: dp
  use cylinder_integrals, only : ljwall_tensor
  implicit none
  real(dp), intent(in) :: k, radiusA, densityA, epsilonratio, sigmaratio
  real(dp) :: local_tensor(3, 3) 
  !! coefficients of r^(-n) for chemical shift, and shielding anisotropy fits:
  !!    n,  chemical shift,       shielding anisotropy
  real(dp), parameter :: coefficients(3, 13) = reshape([ &
       11.0, 2.9972786367663074e9, 5.446624260279856e9, &
       13.0, -5.077950753452642e11, -8.944758704183774e11, &
       15.0, 3.4377090204332965e13, 5.8921865257790164e13, &
       17.0, -1.237329745938399e15, -2.079014276092535e15, &
       19.0, 2.7900236013972024e16, 4.605602208237476e16, &
       21.0, -4.278702675629361e17, -6.938226073887617e17, &
       23.0, 4.646709403875674e18, 7.394035788483907e18, &
       25.0, -3.6303103843431305e19, -5.6603292927874015e19, &
       27.0, 2.0343957863005916e20, 3.10291282943319e20, &
       29.0, -7.997110831817982e20, -1.1910357820108455e21, &
       31.0, 2.0968493616095969e21, 3.0436525190781207e21, &
       33.0, -3.2966096994452803e21, -4.6544214466409094e21, &
       35.0, 2.3521919165514575e21, 3.2234688478153125e21 &
       ], [3, 13])

  local_tensor = ljwall_tensor(k, radiusA, densityA, epsilonratio, &
       sigmaratio, powers = int(coefficients(1, :) + 0.5_dp), &
       isotropic_coeffs = -coefficients(2, :), &
       anisotropy_coeffs = coefficients(3, :))
end function




