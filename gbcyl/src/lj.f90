module lj
use nrtype
use class_parameterizer
use class_parameter_writer
implicit none

interface lj_potential
  module procedure lj1, lj3
end interface

interface lj_init
  module procedure lj_init_parameters, lj_init_wt_reader
end interface

real(dp), save :: sigma_0_ = 1._dp
real(dp), save :: epsilon_0_ = 1._dp

contains

  subroutine lj_init_parameters(epsilon_0, sigma_0)
    real(dp), intent(in) :: epsilon_0, sigma_0
    epsilon_0_ = epsilon_0
    sigma_0_ = sigma_0
  end subroutine

  subroutine lj_init_wt_reader(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'lj_sigma_0', sigma_0_)
    call getparameter(reader, 'lj_epsilon_0', epsilon_0_)
  end subroutine

  subroutine lj_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writeparameter(writer, 'lj_sigma_0', sigma_0_)
    call writeparameter(writer, 'lj_epsilon_0', epsilon_0_)
  end subroutine

  pure function lj1(rij)
    real(dp) :: lj1
    real(dp), intent(in) :: rij
    lj1 = lj3(rij, epsilon_0_, sigma_0_)
  end function

  pure function lj3(rij, epsilon, sigma)
    real(dp) :: lj3
    real(dp), intent(in) :: rij
    real(dp), intent(in) :: epsilon
    real(dp), intent(in) :: sigma
    lj3 = (sigma / rij)**6
    lj3 = lj3 * (lj3 - 1._dp)
    lj3 = 4._dp * epsilon * lj3
  end function

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
