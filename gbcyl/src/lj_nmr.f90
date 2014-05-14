module lj_nmr
!! A helper module for calculating NMR tensors for a spherical particle. 
!!
!! Everything in this module works in the reduced unit system.
!!
use nrtype
use m_constants
use gblj, only: gblj_r, gblj_init
use m_rank2_tensor
use class_poly_box
use particle
use utils
use class_parameterizer
use lj
use particlewall
implicit none


interface 
  pure function ljlj_local(rij) result(t)
    use nrtype, only: dp
    real(dp), intent(in) :: rij
    real(dp) :: t(3, 3)
  end function
end interface

interface 
  pure function gblj_local(x, z) result(t)
    use nrtype, only: dp
    real(dp), intent(in) :: x, z
    real(dp) :: t(3, 3)
  end function
end interface

interface 
  pure function ljwall_local(r, radius) result(t)
    use nrtype, only: dp
    real(dp), intent(in) :: r, radius
    real(dp) :: t(3, 3)
  end function
end interface

contains

subroutine lj_nmr_init(reader)
  type(parameterizer), intent(in) :: reader
  call initptwall(reader)
  call gblj_init(reader)
  call lj_init(reader)
end subroutine


!! Computes the Xe-GB tensor in the local coordinate system in which the z-axis
!! is parallel to the long axis of the GB particle.
!!
!! @p simbox the simulation box in which the particles are located
!! @p gb_particle the Gay-Berne particle
!! @p xe_particle the 131Xe atom
!!
pure function gblj_tensor(tensor, simbox, gb_particle, xe_particle) &
  result(local_tensor)
  interface 
    pure function tensor(x, z)
      use nrtype
      real(dp), intent(in) :: x, z
      real(dp) :: tensor(3, 3)
    end function
  end interface
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: gb_particle
  type(particledat), intent(in) :: xe_particle
  type(rank2_tensor) :: local_tensor
  real(dp) :: rij(3), ui(3)

  rij = minimage(simbox, position(xe_particle) - position(gb_particle))
  ui = orientation(gb_particle)
  local_tensor%axes = gblj_axes(rij, ui)
  local_tensor%components = tensor(dot_product(rij, &
  local_tensor%axes(1:3, 1)), dot_product(rij, local_tensor%axes(1:3, 3)))
end function


!! Returns the local coordinate axis for a GB-LJ pair. 
!!
!! @p rij the vector from GB center to LJ center.
!! @p z the orientation vector of the GB particle.
!! 
!! @return the local coordinate axes. 
!!
pure function gblj_axes(rij, ui) result(axes)
  real(dp), intent(in) :: rij(3), ui(3)
  real(dp) :: axes(3, 3)
  real(dp) :: x(3), y(3), z(3)
  z = ui
  if (dot_product(rij, z) < 0._dp) then
    !! Restrict angle between rij and z below 90 degrees by reverting the z 
    !! axis.
    z = -z
  end if
  !! take the component of rij perpendicular to local z and normalize
  x = rij - dot_product(rij, z) * z
  x = x / sqrt(dot_product(x, x))
  y = cross_product(z, x)
  axes(1:3, 1) = x
  axes(1:3, 2) = y
  axes(1:3, 3) = z
end function


!! Computes the Xe-Xe quadrupolar coupling tensor in the local coordinate
!! system in which the z axis points along the interatomic vector. Unit of
!! the tensor elements is MHz! 
!! 
!! @p simbox the simulation box in which the distances are calculated
!! @p xe the 131Xe particle
!! @p another the other Xe particle
!! 
pure function ljlj_tensor(tensor, simbox, xe, another) result(t)
  interface 
    pure function tensor(r)
      use nrtype
      real(dp), intent(in) :: r
      real(dp) :: tensor(3, 3)
    end function
  end interface
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: xe, another
  real(dp) :: rij(3)
  type(rank2_tensor) :: t
 
  rij = minimage(simbox, position(another) - position(xe))
  t%axes = ljlj_axes(rij)
  t%components = tensor(sqrt(dot_product(rij, rij)))
end function


!! Returns a local axis system for a LJ-LJ pair.
!! 
!! @p rij the center-to-center vector from particle i to particle j.
!!
!! @return the local coordinate axes where z is along rij.
!!
pure function ljlj_axes(rij) result(axes)
  real(dp), intent(in) :: rij(3)
  real(dp) :: axes(3, 3)
  real(dp) :: x(3), y(3), z(3)
  !! Create the local coordinate axes from intermolecular vector rij.
  z = rij / sqrt(dot_product(rij, rij))
  x = [1._dp, 0._dp, 0._dp]
  x = x - dot_product(x, z) * x
  x = x / sqrt(dot_product(x, x))
  y = cross_product(z, x)
  axes(1:3, 1) = x
  axes(1:3, 2) = y
  axes(1:3, 3) = z
end function


!! Computes the Xe-cylindrical cavity quadrupolar coupling / shielding tensor
!! in the local coordinate system in which the z-axis points along the cylinder
!! axis and  x-axis points to the radial direction. Unit of the tensor elements
!! is MHz! The wall consists of smoothly and evenly distributed Lennard-Jones 
!! particles.
!! 
!! @p simbox the simulation box.
!! @p xe the xenon atom.
!!
!! @return the quadrupole coupling / shielding tensor for the 129/131Xe atom
!! and the wall.
!! 
pure function ljwall_tensor(tensor, simbox, xe) result(t)
  interface
    pure function tensor(r, cyl_radius)
      use nrtype
      real(dp), intent(in) :: r, cyl_radius
      real(dp) :: tensor(3, 3)
    end function
  end interface
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: xe
  type(rank2_tensor) :: t
  real(dp) :: r(3), cyl_radius
  r = position(xe)
  cyl_radius = getx(simbox) / 2._dp
  t%components = tensor(sqrt(r(1)**2 + r(2)**2), cyl_radius)
  t%axes = ljwall_axes(r)
end function


!! Returns the local axis system where x-axis points along the radial direction 
!! and goes through the lj particle center. z-axis points along cylinder axis.
!!
!! @p r the position vector of the LJ particle.
!!
!! @return the local coordinate axes.
!!
pure function ljwall_axes(r) result(axes)
  real(dp), intent(in) :: r(3)
  real(dp) :: axes(3, 3)
  real(dp) :: x(3), y(3), z(3)
  !! Create the local coordinate axes.
  z = [0._dp, 0._dp, 1._dp]
  x = r
  x(3) = 0._dp
  x = x / sqrt(dot_product(x, x))
  y = cross_product(z, x)
  axes(1:3, 1) = x
  axes(1:3, 2) = y
  axes(1:3, 3) = z
end function

end module
