module m_shielding
use particle
use class_poly_box
use nrtype
use utils, only: rotate_tensor, cross_product
use gayberne, only: gblj_r, get_gblj_sigma_0
implicit none
external horner
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
!! The shielding is cut off (set to zero) when center-to-center distance
!! exceeds the cutoff value. This is because the Xe-LC fits are reasonable
!! only when r < 5.5.
real(dp), parameter :: cutoff = 5.5_dp 

contains

subroutine avg_shielding(simbox, particles, axes, tensor)
  !! 
  !> Calculates the average Xe shielding in a snapshot as defined in
  !! J. Lintuvuori, M. Straka, and J. Vaara. Nuclear magnetic resonance 
  !! parameters of atomic xenon dissolved in gay-berne model liquid crystal.
  !! Physical Review E, 75(3):031707, MAR 2007.
  !!
  !! @p simbox the simulation volume defining boundary conditions (e.g. 
  !!    perioidic boundaries)
  !! @p particles the array of particles to which the Xe shielding is 
  !!    calculated.
  !! @p axes the coordinate system where the shielding tensor is calculated. 
  !!    axes(:, 1) = x-axis, axes(:, 2) = y-axis, axes(:, 3) = z-axis.
  !! @p tensor the averaged Xe shielding tensor..
  !!
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: particles(:)
  real(dp), intent(in) :: axes(3, 3)
  real(dp), intent(out) :: tensor(3, 3)

  real(dp) :: local_x(3), local_y(3), local_z(3)
  real(dp) :: unit_tensor(3, 3)
  real(dp) :: local_tensor(3, 3)
  real(dp) :: rotated_tensor(3, 3)
  real(dp) :: rij(3)
  integer :: i, j

  tensor = 0._dp
  do i = 1, size(particles)
    if (.not. particles(i)%rod) then
      do j = 1, size(particles)
        if (i == j) cycle
        rij = minimage(simbox, position(particles(j)) - position(particles(i)))
        if(dot_product(rij, rij) > cutoff**2) cycle 
        local_tensor = 0._dp
        local_x = (/1._dp, 0._dp, 0._dp/)
        local_y = (/0._dp, 1._dp, 0._dp/)
        local_z = (/0._dp, 0._dp, 1._dp/)
        unit_tensor = reshape((/local_x, local_y, local_z/), (/3, 3/))
        !! Calculate tensor components in a coordinate frame fixed to the molecule pair j, i.
        if (particles(j)%rod) then 
          call gblj_shielding(particles(j), rij, local_tensor)
          !! determine rotation angles from particles(j)%(ux, uy, uz) and rij
          local_z = orientation(particles(j))
          !! take the component of rij perpendicular to local z and normalize
          local_x = rij - dot_product(rij, local_z) * local_z
          local_x = local_x / sqrt(dot_product(local_x, local_x))
        else 
          call ljlj_shielding(rij, local_tensor)
          local_z = rij / sqrt(dot_product(rij, rij))
          local_x = local_x - dot_product(local_x, local_z) * local_x
          local_x = local_x / sqrt(dot_product(local_x, local_x))
          !! determine rotation angles from rij
        end if
        local_y = cross_product(local_z, local_x)
        !write(*, *) local_y, dot_product(local_y, local_y)
        !! rotate local_tensor to laboratory coordinates.
        ! if laboratory coordinate frame is the frame of the simbox then new_axes = unit_tensor
        call rotate_tensor(local_tensor, old_axes = &
          reshape((/local_x, local_y, local_z/), (/3, 3/)), &
          new_axes = axes, rotated = rotated_tensor)
        !write(*, *) "Tr(local) = ", local_tensor(1, 1) + local_tensor(2, 2) + local_tensor(3, 3)
        !write(*, *) "Tr(rotated) = ", rotated_tensor(1, 1) + rotated_tensor(2, 2) + rotated_tensor(3, 3)
  
        !! Add contribution of particle j to the total shielding
        forall(i = 1:3, j = 1:3) tensor(i, j) = tensor(i, j) + rotated_tensor(i, j)      
      end do 
    end if
  end do
  !! Divide by number of lj particles
  tensor = tensor / real(count(.not. particles%rod), dp)
end subroutine


!> Calculates the NMR shielding tensor for a GB-Xe pair. The shielding tensor
!! is calculated in a local coordinate system in which the orientation vector
!! of the GB particle defines the z-axis and the perpendicular component of the
!! interparticle distance vector defines the x-axis.
!! 
!! @p gb_particle the Gay-Berne particle with respect to which the shielding is
!! calculated. 
!! @p lj_particle the Xe atom to which the shielding is calculated. 
!! @p rij the center to center (minimum image) distance vector of the particles
!! @p local_tensor the shielding tensor calculated in a coordinate system 
!! defined by the geometry of the particle pair.
!! 
subroutine gblj_shielding(gb_particle, rij, local_tensor)
  type(particledat), intent(in) :: gb_particle
  real(dp), intent(in) :: rij(3)
  real(dp), intent(out) :: local_tensor(3, 3)
  real(dp) :: co, si
  real(dp) :: urij(3)
  real(dp) :: r
  urij = rij / sqrt(dot_product(rij, rij))
  !! cosine and sine of the angle between the orientation of the gb_particle
  !! and the interparticle vector rij:
  co = dot_product(orientation(gb_particle), urij)
  si = sqrt(1._dp - co**2) 
  !! We can always say that sine is positive because the angle between is
  !! smaller or equal to pi. For definition of the angle see the paper by
  !! Lintuvuori et al., Figure 1.
  r = gblj_r(orientation(gb_particle), rij) * get_gblj_sigma_0()
  local_tensor = 0._dp
  local_tensor(1, 1) = sigma(r, s_xx) * si**2 + &
    sigma(r, e_perpendicular) * co**2 !! sigma_xx
  local_tensor(2, 2) = sigma(r, e_perpendicular) !! sigma_yy
  local_tensor(3, 3) = sigma(r, s_parallel) * si**2 !! sigma_zz
  local_tensor(1, 3) = sigma(r, xz) * si * co !! sigma_xz
  local_tensor(3, 1) = sigma(r, zx) * si * co !! sigma_zx
end subroutine 


!> Calculates the NMR shielding tensor for two Xe atoms described as LJ 
!! particles.
!!
!! @p lj_1 the LJ particle with respect to which the shielding is calculated.
!! @p lj_2 the LJ particle to which the shielding is calculated.
!! @p rij the center to center (minimum image) distance vector of the particles
!! @p local_tensor the tensor calculated in coordinates where rij defines the 
!! z-axis.
!!
subroutine ljlj_shielding(rij, local_tensor)
  real(dp), intent(in) :: rij(3)
  real(dp), intent(out) :: local_tensor(3, 3)
  real(dp) :: r_aengstroms
  local_tensor = 0._dp
  r_aengstroms = sqrt(dot_product(rij, rij)) * 4.5_dp !! 1 GB unit is here 4.5 Å 
  if (r_aengstroms < cutoff) then 
    local_tensor(1, 1) = sigma(r_aengstroms, xexe_isotropic) - sigma(r_aengstroms, anisotropy_xexe) / 3._dp 
    local_tensor(2, 2) = sigma(r_aengstroms, xexe_isotropic) - sigma(r_aengstroms, anisotropy_xexe) / 3._dp 
    local_tensor(3, 3) = sigma(r_aengstroms, xexe_isotropic) + 2._dp/3._dp * sigma(r_aengstroms, anisotropy_xexe)
  end if
end subroutine 


!> Implements the function (17) from the paper by Lintuvuori et al. 
!! 
!! @p r the distance
!! @p a the parameter A and polynomial coefficients in an array [A, p_A0, p_A1,...]
!! @p b the parameter B and polynomial coefficients in an array [B, p_B0, p_B1,...] 
!!
function sigma(r, a, b) result(res)
  real(dp), intent(in) :: R
  real(dp), intent(in) :: a(:)
  real(dp), intent(in), optional :: b(:)
  real(dp) :: res
  interface
  pure function horner(a, n, x)
    use nrtype
    real(dp), dimension(:), intent(in) :: a
    integer, intent(in) :: n 
    real(dp), intent(in) :: x
    real(dp) :: horner
  end function horner
  end interface
  if (r > cutoff) write(*, *) 'r=', r, "WTF?!!"
  res = a(1) / (r**horner(a(size(a):2:-1), size(a) - 1, r))
  !write(*, *) horner(a(size(a):2:-1), size(a) - 1, r)
  if (present(b)) res = res + b(1) / r**horner(b(size(b):2:-1), size(b) - 1, R)
end function

end module
