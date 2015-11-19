!! Computes the quadrupole coupling tensor for a 131Xe atom
!! inside a cylindrical cavity with radius @p radiusA. Cavity wall is
!! infinitely thick and consists of smoothly and evenly distributed
!! Lennard-Jones particles. The local coordinate system has z as the
!! cylinder axis and x pointing to the radial direction. y-axis is
!! perpendicular to both. The @p epsilonratio is used to scale the strength
!! of the shielding with respect to the Xe-Xe shielding. @p sigmaratio is
!! used to scale the distance to the / curvature of the wall compared to a
!! wall consisting of Xe atoms. 
!! 
!! @p k the ratio (distance of Xe from cylinder axis) / 
!!    (inner radius of the cylindrical cavity)
!! @p radiusA the inner radius of the cavity in aengstroms
!! @p densityA the density of the wall in 1/aengstroms^3
!! @p epsilonratio = xewall_epsilon_0 / xexe_epsilon_0
!! @p sigmaratio = xewall_sigma_0 / xexe_sigma_0
!!
!! @return local_tensor the quadrupole coupling tensor (kHz) in the local
!!  coordinate system.
!!
function xewall_qcoupling(k, radiusA, densityA, epsilonratio, &
  sigmaratio) result(local_tensor)
  use num_kind, only: dp
  use cylinder_integrals, only : ljwall_tensor
  real(dp), intent(in) :: k, radiusA, densityA, epsilonratio, sigmaratio
  real(dp) :: local_tensor(3, 3) 
  !! coefficients of r^(-n) for the "parallel" component of the Xe-Xe 
  !! quadrupole coupling tensor. 
  !!    n,  khi_parallel
  real(dp), parameter :: khiparall(2, 13) = reshape([ & 
       11., -7.678400838773115e11, &
       13., 9.558665403725903e13, & 
       15., -5.988138037657408e15, &
       17., 2.3324762697422534e17, &
       19., -6.017812289055779e18, &
       21., 1.069234232617967e20, &
       23., -1.340929792227524e21, &
       25., 1.1996416647616803e22, &
       27., -7.625396244715735e22, &
       29., 3.3690967454200235e23, &
       31., -9.847039785026325e23, &
       33., 1.713090762707716e24, &
       35., -1.3438977967328407e24 &
       ], [2, 13])
  local_tensor = ljwall_tensor(k, radiusA, densityA, epsilonratio, sigmaratio,&
       powers = int(khiparall(1, :) + 0.5_dp), &
       anisotropy_coeffs = 1.5_dp * khiparall(2, :))
end function



