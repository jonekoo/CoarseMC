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
   pure function zhangI(m, k, rc)
     use num_kind, only: dp
     implicit none
     integer, intent(in) :: m
     real(dp), intent(in) :: k, rc
     real(dp) :: zhangI
   end function zhangI
end interface

interface
   pure function d_zhangI(m, k, rc)
     use num_kind, only: dp
     implicit none
     integer, intent(in) :: m
     real(dp), intent(in) :: k, rc
     real(dp) :: d_zhangI
   end function d_zhangI
end interface

interface 
   real(c_double) pure function hyp2f1(a, b, c, z) bind(C, name='hyp2f1')
     use, intrinsic :: iso_c_binding
     real(c_double), intent(in), value :: a, b, c, z
   end function hyp2f1
end interface

interface
   real(c_double) pure function gamma(x) bind(C, name='gamma')
     use, intrinsic :: iso_c_binding
     real(c_double), intent(in), value :: x
   end function gamma
end interface

end module

