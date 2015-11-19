!> Implements the interaction for a spherical particle inside a
!! cylindrical cavity. The particle interacts with the elements in the
!! cavity wall through the Lennard-Jones 12-6 potential. The total
!! interaction is integrated over the volume of the infinitely thick
!! wall of an infinitely long cavity.
!!
!! @see X. Zhang and W. Wang and G. Jiang, "A potential model for
!! interaction between the Lennard-Jones cylindrical wall and fluid
!! molecules". Fluid Phase Equilibria, 218(2), 239-246, 2004.


module ljcylinder
use num_kind
implicit none

interface 
   !> Returns the repulsive term of the interaction between a LJ
   !! particle and the smooth wall of the long cylindrical cavity.
   !! @p k = r / rc, where r is the distance from the cavity center
   !! and @p rc is the inner radius of the cavity. The routine works in
   !! reduced LJ units, where distances are in sigma_0:s and energies
   !! in epsilon_0:s.
   pure function repwall2(k, rc)
     use num_kind
     real(dp), intent(in) :: k, rc
     real(dp) :: repwall2
   end function repwall2

   !> Returns the attractive term of the interaction between a LJ
   !! particle and the smooth wall of the long cylindrical cavity.
   !! @p k = r / rc, where r is the distance from the cavity center
   !! and @p rc is the inner radius of the cavity. The routine works in
   !! reduced LJ units, where distances are in sigma_0:s and energies
   !! in epsilon_0:s.
   pure function attwall2(k, rc)
     use num_kind
     real(dp), intent(in) :: k, rc
     real(dp) :: attwall2
   end function attwall2

   !> Returns the repulsive term of the force acting on a LJ
   !! particle by the smooth wall of the long cylindrical cavity.
   !! @p k = r / rc, where r is the distance from the cavity center
   !! and @p rc is the inner radius of the cavity. The routine works in
   !! reduced LJ units, where distances are in sigma_0:s and energies
   !! in epsilon_0:s.
   pure real(dp) function d_repwall(k, rc)
     use num_kind
     real(dp), intent(in) :: k, rc
   end function d_repwall

   !> Returns the attractive term of the force acting on a LJ
   !! particle by the smooth wall of the long cylindrical cavity.
   !! @p k = r / rc, where r is the distance from the cavity center
   !! and @p rc is the inner radius of the cavity. The routine works in
   !! reduced LJ units, where distances are in sigma_0:s and energies
   !! in epsilon_0:s.
   pure real(dp) function d_attwall(k, rc)
     use num_kind
     real(dp), intent(in) :: k, rc
   end function d_attwall
end interface

contains

!> Returns the potential energy for a Lennard-Jones (LJ) site with
!! respect to a cylindrical wall consisting of smoothly and evenly
!! distributed LJ particles.
!! 
!! @param eps the parameter for LJ-particle - wall well-depth.
!! @param density the wall density.
!! @param sigma the parameter setting the contact distance between a
!!        wall element and the LJ site.
!! @param alpha controls the amount of attraction with respect to
!!        repulsion. alpha=1 is the regular LJ interaction while
!!        alpha=0 makes the interaction purely repulsive.
!! @param r the distance of LJ site center from the cylinder axis.
!! @param rwall the radius of the cylinder.
!! 
pure function ljcylinderpotential(eps, density, sigma, alpha, r, rwall) &
     result(potential)
  real(dp), intent(in) :: eps, density, sigma, alpha, r, rwall
  real(dp) :: potential
  potential = eps * density * sigma**3 * (repwall2(r / rwall, rwall / sigma) &
       - alpha * attwall2(r / rwall, rwall / sigma))
end function

!> Returns the value of force acting on a LJ particle due to a
!! cylindrical wall consisting of smoothly and evenly distributed LJ
!! particles. The direction of the force is from the center of the
!! cavity towards the wall. 
!!
!! @param eps the parameter for LJ-particle - wall well-depth
!! @param density the wall density
!! @param sigma the parameter setting the contact distance for the
!!        wall element and the LJ site.
!! @param alpha controls the amount of attraction with respect to
!!        repulsion. alpha=1 is the regular LJ interaction while
!!        alpha=0 makes the interaction purely repulsive.
!! @param r the distance of LJ site center from the cylinder axis.
!! @param rwall the radius of the cylinder.
!! 
pure function ljcylinderforce(eps, density, sigma, alpha, r, rwall)
  real(dp), intent(in) :: eps, density, sigma, alpha, r, rwall 
  real(dp) :: ljcylinderforce
  ljcylinderforce = eps * density * sigma**3 * &
       (d_repwall(r / rwall, rwall / sigma) - &
       alpha * d_attwall(r / rwall, rwall / sigma)) / rwall
end function

end module
