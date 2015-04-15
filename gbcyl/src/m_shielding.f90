module m_shielding
use nrtype
use utils, only: rotate_tensor, cross_product, horner
use m_gblj
use m_constants
use class_parameterizer !! for initialization
use m_lj
use particlewall
implicit none
!!
!! This module implements the calculation of Xe NMR shielding for a 
!! configuration or "snapshot" of a mixture of Xe and GB particles as defined
!! in J. Lintuvuori, M. Straka, and J. Vaara. Nuclear magnetic resonance 
!! parameters of atomic xenon dissolved in gay-berne model liquid crystal.
!! Physical Review E, 75(3):031707, MAR 2007. The module uses the newer and 
!! better Xe-Xe shielding parameters from M. Hanni, P. Lantto, M. Iliaš,
!! H. Jørgen Jensen, and J. Vaara. Relativistic effects in the intermolecular
!! interaction-induced nuclear magnetic resonance parameters of xenon dimer. 
!! The Journal of Chemical Physics, 127(16):164313, 2007.
!!
!! The functional form is the same for all shielding contributions and
!! implemented as the function sigma in this module. The parameters given to it
!! are listed in the tables/arrays below.
!!
!!
!! Xe-LC NMR shielding fit parameters. (Excerpt from) Table II from the paper 
!! by Lintuvuori et al. Reduced Gay-Berne units.
!!
!!                    A (ppm)            p_A,0               p_A,1               p_A,2              p_A,3
real(dp), parameter :: &
s_xx(4) =            [-1.2026297065e1_dp, 3.3134189498_dp,    9.0507886417_dp,   -1.6025243568_dp],& 
e_perpendicular(4) = [-1.1672648797e1_dp, 3.9997702047e-2_dp, 1.1636011366e1_dp, -1.5288591743_dp],&
s_parallel(5) =      [-5.2671604883e1_dp, 4.5165162489_dp,    1.5171022733e1_dp, -2.3590723634e1_dp, 9.9306821259_dp],&
xz(3) =              [ 3.1228263091e1_dp, 5.3141730331_dp,    6.6334674074e-1_dp],&    
zx(3) =              [ 3.3141615175e1_dp, 4.8095621340_dp,    7.0308516233e-1_dp]
!!
!! The Xe-LC fits are reasonable only when r < 5.5. When center-to-center 
!! distance exceeds this value, the shielding is set to zero.
real(dp), parameter :: gbxe_cutoff = 5.5_dp 
!!
!!
!! Fit parameters for the Xe-Xe shielding from the paper by M. Hanni, et al. 
!! Xe-Xe Distance has to be in aengstroms (Å) when using these! The function 
!! ljlj_shielding contains a hardcoded definition 1 GB length unit = 4.5 Å. 
!!
!!                    A (Å ppm)        p_0                 p_1 (Å^-1)        p_2 (Å^-2)        p_3 (Å^-3)
real(dp), parameter :: &
xexe_isotropic(3) =  [-3329.92155454_dp, -1.73732128_dp,      1.11440532_dp],&                   !! (*)
anisotropy_xexe(4) = [ 6274.20005258_dp, -1.20014259_dp,      0.93600894_dp,      0.02009022_dp] !! (**)
!!
!! (*) The isotropic part is here actually the fit of the chemical shift. The 
!! parameter A for the chemical shift just does not have the minus sign:
!! chemical shift = reference shielding - shielding.
!! When the reference is zero, we get the shielding = -chemical shift.
!!
!! (**) This is a fit for anisotropy of shielding, not the anisotropic part of
!! shielding.
!!
!! The quantum chemical calculations of Hanni were made in the region
!! 3...13.5 Å. When the center-to-center distance of two Xe atoms exceeds 
!! 13.5 Å the shielding is set to zero.
real(dp), parameter :: xexe_cutoff_A = 13.5_dp

type(gblj_potential) :: gblj
type(lj_potential) :: lj

interface xewall_shielding
  pure function xewall_shielding(k, radiusA, densityA, epsilonratio, &
    sigmaratio) result(local_tensor)
    use nrtype
    real(dp), intent(in) :: k, radiusA, densityA, epsilonratio, sigmaratio
    real(dp) :: local_tensor(3, 3)
  end function
end interface

contains


!> Initializes the module and the modules this module depends on.
!! 
!! @param reader the object responsible for getting the parameters from e.g.
!! an input file.
!!
subroutine init_shielding(reader)
  type(parameterizer), intent(inout) :: reader
  !! The modules below are needed for the ljwall_shielding calculation
  call particlewall_init(reader)
  lj = lj_potential(reader)   
  gblj = gblj_potential(reader)
end subroutine



!> Calculates the NMR shielding tensor for a GB-Xe pair. The shielding tensor
!! is calculated in a local coordinate system in which the orientation vector
!! of the GB particle defines the z-axis and the perpendicular component of the
!! interparticle distance vector defines the x-axis.
!! 
!! @param x,z the x and z coordinates of the Xe atom in the local axis system.
!!
!! @return the contribution of one GB particle to the 129Xe shielding tensor.
!!
pure function gbxe_shielding_local(x, z) result(local_tensor)
  real(dp), intent(in) :: x, z
  real(dp) :: local_tensor(3, 3)
  real(dp) :: co, si
  real(dp) :: r
  local_tensor = 0._dp
  if (z**2 + x**2 < gbxe_cutoff**2) then
     !! cosine and sine of the angle between the orientation of the gb_particle
     !! and the interparticle vector rij:
     co = z / sqrt(z**2 + x**2)
     si = sqrt(1._dp - co**2)
     r = gblj%r([0._dp, 0._dp, 1._dp], [x, 0._dp, z])
     local_tensor = 0._dp
     local_tensor(1, 1) = sigma(r, s_xx) * si**2 + &
          sigma(r, e_perpendicular) * co**2 !! sigma_xx
     local_tensor(2, 2) = sigma(r, e_perpendicular) !! sigma_yy
     local_tensor(3, 3) = sigma(r, s_parallel) * si**2 !! sigma_zz
     local_tensor(1, 3) = sigma(r, xz) * si * co !! sigma_xz
     local_tensor(3, 1) = sigma(r, zx) * si * co !! sigma_zx
  end if
end function


!> Calculates the NMR shielding tensor for two 129/131Xe atoms. 
!!
!! @param rij is the center-to-center distance of the atoms.
!!
pure function xexe_shielding_local(rij) result(local_tensor)
  real(dp), intent(in) :: rij
  real(dp) :: local_tensor(3, 3)
  real(dp) :: r_aengstroms
  local_tensor = 0._dp
  r_aengstroms = rij * sigma0_aengstroms
  if (r_aengstroms < xexe_cutoff_A) then 
     local_tensor(1, 1) = sigma(r_aengstroms, xexe_isotropic) - &
          sigma(r_aengstroms, anisotropy_xexe) / 3._dp 
     local_tensor(2, 2) = sigma(r_aengstroms, xexe_isotropic) - &
          sigma(r_aengstroms, anisotropy_xexe) / 3._dp 
     local_tensor(3, 3) = sigma(r_aengstroms, xexe_isotropic) + 2._dp/3._dp * &
          sigma(r_aengstroms, anisotropy_xexe)
  end if
end function


!> Implements the function (17) from the paper by Lintuvuori et al. 
!! 
!! @param r the distance
!! @param a the parameter A and polynomial coefficients in an array 
!!    [A, p_A0, p_A1,...]
!! @param b the parameter B and polynomial coefficients in an array 
!!    [B, p_B0, p_B1,...] 
!!
pure function sigma(r, a, b) result(res)
  real(dp), intent(in) :: R
  real(dp), intent(in) :: a(:)
  real(dp), intent(in), optional :: b(:)
  real(dp) :: res
  res = a(1) / (r**horner(a(size(a):2:-1), r))
  if (present(b)) res = res + b(1) / r**horner(b(size(b):2:-1), R)
end function


!> Returns the contribution of a smooth, cylindrical Lennard-Jones wall to the
!! nuclear shielding of a Xe atom. Transforms the units to aengstroms (Å) 
!! before passing the actual computation to the function xewall_shielding. 
!!
!! The computation is performed in the local coordinate system in which the 
!! z-axis is parallel to the cylinder axis and x axis points to the radial 
!! direction.
!!
!! @see module particlewall
!!
!! @param r the distance of the Xe atom from the cylinder axis.
!! @param radius the inner radius of the cylindrical wall.
!!
!! @return the Xe-wall shielding tensor.
!!
pure function xewall_shielding_local(r, radius) result(local_tensor)
  use particlewall, only: ljwall_sigma_0 => sigwall_lj, ljwall_epsilon_0 => &
  epswall_lj, wall_density
  real(dp), intent(in) :: radius
  real(dp), intent(in) :: r
  real(dp) :: local_tensor(3, 3) 
  local_tensor = 0._dp
  local_tensor = xewall_shielding(r / radius, &
    radius * sigma0_aengstroms, &
    wall_density / sigma0_aengstroms**3, ljwall_epsilon_0 / lj%epsilon_0, &
    ljwall_sigma_0 / lj%sigma_0)
end function

end module
