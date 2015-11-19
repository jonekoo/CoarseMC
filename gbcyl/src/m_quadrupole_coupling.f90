module m_quadrupole_coupling
use num_kind
use m_constants
use m_gblj
use m_shielding, only: sigma
use utils
use class_parameterizer
use m_lj
use particlewall
implicit none

!! Lintuvuori et al. Phys. Rev. E 75, 031707 (2007), excerpt from TABLE II: 
!! Parameters of 131Xe quadrupole coupling tensors when interacting with the
!! model liquid crystal. The unit of A and B is MHz. The $p_i$ are in units of
!! $\sigma_0^{-i}$.  
!!
!!                A                 p_A,0            p_A,1           p_A,2            p_A,3
real(dp), parameter :: &
xx_s(4)       = [ 2.9222586183,     3.2835741349,    3.6896880124,  -6.4687585514e-1], &                 !! $\chi_{xx}^s$ 
a_parall_e(4) = [-5.2004492490e-1,  1.6846133021e1, -1.0664334268e1, 2.3313944890], &                    !! $\chi_{||}^e$ 
parall_s(5)   = [-3.8213603446e-1,  9.9725303666,   -5.2124031847,  -6.8656465241e-1, 6.8925841175e-1],& !! $\chi_{||}^s$
a_xz(4)       = [-7.6955322406e-6, -3.8558980871e1, -3.4192991767e1, 2.8627512246e1]                     !! $\chi_{xz}$

!!                B                 p_B,0            p_B,1           p_B,2           
real(dp), parameter :: & 
b_parall_e(4) = [ 2.9010256384e-1,  2.6871761243e1, -4.6117061897e1, 3.7140396050e1], &                  !! $\chi_{||}^e$
b_xz(4)       = [ 1.3445614681,     6.5490556256e1, -1.7162263010e2, 1.2463647343e2]                     !! $\chi_{xz}$

!! The coupling is cut off (set to zero) when center-to-center distance
!! exceeds the cutoff value. This is because the Xe-LC fits are reasonable
!! only when r < 5.5.
real(dp), parameter :: gbxe_cutoff = 5.5_dp 



!! Hanni et al. J. Chem. Phys. 127, 164313 (2007), excerpt from TABLE I: 
!! Fitting parameters for 131Xe nuclear quadrupole coupling $\chi_{||}(R)$,
!! in kHz(!), in Xe2.
!!
!!                   A/B                    p_0/q_0       p_1/q_1      p_2/q_2
real(dp), parameter :: &
a_xexe_parall(4) = [ 1620336.53452898,      2.12902635,  -0.13734912,  0.16125823],& !! $\chi_{||}$ (A-term)
b_xexe_parall(4) = [ 13488314876.57440000,  16.12455833, -1.45189255,  0.08429018]   !! $\chi_{||}$ (B-term)
!!
!! 
!! The quantum chemical calculations of Hanni were made in the region
!! r = 3...13.5 Å. When the center-to-center distance of two Xe atoms exceeds
!! 13.5 Å, the coupling is set to zero.
real(dp), parameter :: xexe_cutoff_A = 13.5_dp

type(gblj_potential) :: gblj
type(lj_potential) :: lj

interface 
  pure function xewall_qcoupling(k, radiusA, densityA, epsilonratio, &
    sigmaratio) result(local_tensor)
    use num_kind, only: dp
    real(dp), intent(in) :: k, radiusA, densityA, epsilonratio, sigmaratio
    real(dp) :: local_tensor(3, 3)
  end function
end interface

contains

subroutine qcoup_init(reader)
  type(parameterizer), intent(in) :: reader
  call particlewall_init(reader)
  !call gblj_init(reader)
  gblj = gblj_potential(reader)
  lj = lj_potential(reader)
end subroutine


!! Calculates the components of the quadrupole coupling tensor components in 
!! the local coordinate system, where z-axis is along the long axis of the 
!! GB-ellipsoid and x-axis is the projection of the vector rij connecting the 
!! centers of the particles, perpendicular to z. Direction is chosen so that 
!! the angle between rij and z is always larger than zero and less than pi/2.
!! For illustration see the paper by Lintuvuori et al.
!!
!! @p x is the x-coordinate of the 131Xe atom, x >= 0. 
!! @p z is the z-coordinate of the 131Xe atom, z >= 0. 
!!
pure function gbxe_coupling_local(x, z) result(coupling)
  real(dp), intent(in) :: x, z  
  real(dp) :: coupling(3, 3)
  real(dp) :: r
  real(dp) :: co, si
  real(dp) :: xx, zz, xz

  if (x**2 + z**2 > gbxe_cutoff**2) then
    !! cutoff is needed because at least chi_xx_s explodes near r = 10.
    coupling = 0._dp
    return
  end if

  !! cosine and sine of the angle between the orientation of the gb_particle
  !! and the interparticle vector rij:
  co = z / sqrt(z**2 + x**2)
  si = sqrt(1._dp - co**2) 

  r = gblj%r([0._dp, 0._dp, 1._dp], [x, 0._dp, z])

  !! chi tensor coupling:
  xx = sigma(r, xx_s) * si**2 &
       - 0.5_dp * (sigma(r, a_parall_e) + sigma(r, b_parall_e)) * co**2
  zz = sigma(r, parall_s) * si**2 &
       + (sigma(r, a_parall_e) + sigma(r, b_parall_e)) * co**2
  xz = (sigma(r, a_xz) + sigma(r, b_xz)) * si * co
  coupling = 0._dp
  coupling(1, 1) = xx
  coupling(2, 2) = -(xx + zz)
  coupling(3, 3) = zz
  coupling(1, 3) = xz
  coupling(3, 1) = xz
end function


!! Returns the quadrupole coupling tensor for a Xe dimer in coordinates where
!! z-axis points from center to center.
!! 
!! @p r the center-to-center distance of the xenons.
!! 
!! @return the coupling tensor in MHz.
!!
pure function xexe_coupling_local(r) result(t)
  real(dp), intent(in) :: r
  real(dp) :: t(3, 3)
  real(dp) :: r_aengstroms
  t = 0._dp
  !! Convert distance to aengstroms:
  r_aengstroms = r * sigma0_aengstroms
  if (r_aengstroms < xexe_cutoff_A) then
     t(3, 3) = sigma(r_aengstroms, a_xexe_parall) - sigma(r_aengstroms, b_xexe_parall)
     t(1, 1) = -0.5_dp * t(3, 3)
     t(2, 2) = -0.5_dp * t(3, 3)
     t = t * 0.001_dp !! Conversion to MHz
  end if
end function

!! Computes the Xe-cylindrical cavity quadrupolar coupling tensor in the local
!! coordinate system in which the z-axis points along the cylinder axis and
!! x-axis points to the radial direction. Unit of the tensor elements is MHz! 
!! All input must be expressed in reduced units! The wall consists of smoothly
!! and evenly distributed Lennard-Jones particles.
!! 
!! @p r the Xe:s position vector.
!! @p cyl_radius the radius of the cavity.
!! @p rho the density of the cavity wall. 
!! @p local_tensor the quadrupolar coupling tensor in the local coordinates.
!! 
pure function xewall_coupling_local(r, cyl_radius) result(t)
  use particlewall, only: ljwall_epsilon_0 => epswall_lj, &
    ljwall_sigma_0 => sigwall_lj, wall_density
  real(dp), intent(in) :: r, cyl_radius
  real(dp) :: t(3, 3)
  t = 0._dp

  t = xewall_qcoupling(r / cyl_radius, cyl_radius * sigma0_aengstroms, &
    wall_density / sigma0_aengstroms**3, ljwall_epsilon_0 / lj%epsilon_0, &
    ljwall_sigma_0 / lj%sigma_0)
  t = 0.001_dp * t !! Conversion to MHz

end function

end module
