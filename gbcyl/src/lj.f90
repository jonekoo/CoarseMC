!> Implements the Lennard-Jones 12-6 potential and the force derived from it.
module lj
use nrtype
use class_parameterizer
use class_parameter_writer
implicit none


!> A generic interface for computing the LJ12-6 potential value.
interface lj_potential
  module procedure lj1, lj3
end interface

!> Initializes the module.
interface lj_init
  module procedure lj_init_parameters, lj_init_wt_reader
end interface

!> The range parameter which gives the zero-crossing distance of the
!! potential.
real(dp), save :: sigma_0_ = 1._dp

!> The well-depth of the potential.
real(dp), save :: epsilon_0_ = 1._dp

contains

!> Initializes the LJ 12-6 potential with the well-depth @p epsilon_0
!! and contact distance @p sigma_0.
subroutine lj_init_parameters(epsilon_0, sigma_0)
  real(dp), intent(in) :: epsilon_0, sigma_0
  epsilon_0_ = epsilon_0
  sigma_0_ = sigma_0
end subroutine

!> Initialize the module using parameters given by @p reader.
subroutine lj_init_wt_reader(reader)
  type(parameterizer), intent(in) :: reader
  call getparameter(reader, 'lj_sigma_0', sigma_0_)
  call getparameter(reader, 'lj_epsilon_0', epsilon_0_)
end subroutine

!> Write the module parameters using the output unit and format defined
!! by @p writer.
subroutine lj_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call writeparameter(writer, 'lj_sigma_0', sigma_0_)
  call writeparameter(writer, 'lj_epsilon_0', epsilon_0_)
end subroutine

!> Returns the LJ 12-6 potential with the parameterization set with
!! lj_init and a given internuclear vector @p rij.
pure function lj1(rij)
  real(dp) :: lj1
  real(dp), intent(in) :: rij
  lj1 = lj3(rij, epsilon_0_, sigma_0_)
end function lj1

!> Returns the LJ 12-6 potential value with given well-depth @p epsilon, 
!! range parameter @p sigma and the internuclear vector @p rij from
!! particle i to particle j. 
pure function lj3(rij, epsilon, sigma)
  real(dp) :: lj3
  real(dp), intent(in) :: rij
  real(dp), intent(in) :: epsilon
  real(dp), intent(in) :: sigma
  lj3 = (sigma / rij)**6
  lj3 = lj3 * (lj3 - 1._dp)
  lj3 = 4._dp * epsilon * lj3
end function lj3

!> Returns the force exerted on the particle j by particle i when the
!! interparticle vector from i to j is @p rij.
pure function lj_force(rij)
  real(dp), intent(in) :: rij(3)
  real(dp) :: lj_force(3)
  real(dp) :: urij(3)
  real(dp) :: rijabs
  rijabs = sqrt(dot_product(rij, rij))
  urij = rij/rijabs
  lj_force = 24._dp * epsilon_0_ * sigma_0_**6/rijabs**7 * &
    (2._dp * sigma_0_**6/rijabs**6 - 1._dp) * urij
end function

end module
