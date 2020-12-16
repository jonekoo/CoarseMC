!> Interfaces of the integrals used in the particle-wall interactions
!! for a Lennard-Jones particle inside a cylindrical cavity. The walls
!! of the cavity consist of smoothly and evenly distributed LJ
!! particles.
module cylinder_integrals

interface
   function ljwall_tensor(k, radiusA, densityA, epsilonratio, &
        sigmaratio, powers, isotropic_coeffs, anisotropy_coeffs) &
        result(local_tensor)
     use num_kind, only: dp
     implicit none
     real(dp), intent(in) :: k, radiusA, densityA, epsilonratio, sigmaratio
     integer, intent(in) :: powers(:)
     real(dp), intent(in), optional :: isotropic_coeffs(:)
     real(dp), intent(in), optional :: anisotropy_coeffs(:)
     real(dp) :: local_tensor(3, 3) 
   end function ljwall_tensor
end interface

interface
   !> Computes numerically the integral I_m(k, R), k = r/R, in the paper
   !! X. Zhang, W. Wang, and G. Jiang, Fluid Phase Equilibria 2004,
   !! 218(2), 239-246.
   pure function zhangI(m, k, rc)
     use num_kind, only: dp
     implicit none
     integer, intent(in) :: m
     real(dp), intent(in) :: k, rc
     real(dp) :: zhangI
   end function zhangI
end interface

interface
   !> Computes numerically the derivative d I_m(k, R) / d k of the
   !! integral I_m(k, R), k = r/R, in the paper X. Zhang, W. Wang, and
   !! G. Jiang, Fluid Phase Equilibria 2004, 218(2), 239-246.
   pure function d_zhangI(m, k, rc)
     use num_kind, only: dp
     implicit none
     integer, intent(in) :: m
     real(dp), intent(in) :: k, rc
     real(dp) :: d_zhangI
   end function d_zhangI
end interface

interface
   !> Computes numerically the value of the Hypergeometric function
   !! _2F_1(a, b, c, z), as defined in the Handbook of Mathematical
   !! Functions (1972) by Abramovitz and Stegun.
   real(c_double) pure function hyp2f1(a, b, c, z) bind(C, name='hyp2f1')
     use, intrinsic :: iso_c_binding
     real(c_double), intent(in), value :: a, b, c, z
   end function hyp2f1
end interface

interface
   !> Computes the value of the gamma function at @p x.
   real(c_double) pure function gamma(x) bind(C, name='gamma')
     use, intrinsic :: iso_c_binding
     real(c_double), intent(in), value :: x
   end function gamma
end interface

end module

