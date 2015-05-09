!> Implements the Lennard-Jones 12-6 potential and the force derived from it.
module m_lj
  use nrtype
  use class_parameterizer
  use class_parameter_writer
  use m_point_interaction, only: pointpoint_potential
  use m_point, only: point
  implicit none
  
  !> Initializes the module.
  interface lj_init
     module procedure lj_init_parameters, lj_init_wt_reader
  end interface lj_init
  
  type, extends(pointpoint_potential) :: lj_potential
     !> The range parameter which gives the zero-crossing distance of the
     !! potential.
     real(dp) :: sigma_0 = 1.
     
     !> The well-depth of the potential.
     real(dp) :: epsilon_0 = 1.
   contains
     procedure :: lj1
     procedure :: potential => lj_potential_potential
     procedure :: force => lj_force
     procedure :: writeparameters => lj_writeparameters
  end type lj_potential
  
  interface lj_potential
     module procedure lj_init_parameters, lj_init_wt_reader
  end interface lj_potential
  
  
contains

  subroutine lj_potential_potential(this, particle_i, particle_j, res, err)
    class(lj_potential), intent(in) :: this
    type(point), intent(in) :: particle_i, particle_j
    real(dp), intent(out) :: res
    integer, intent(out) :: err
    logical :: overlap
    err = 0
    call this%lj1(norm2(particle_j%position() - particle_i%position()), &
         res, overlap)
    if (overlap) err = -1
  end subroutine lj_potential_potential
  
  !> Initializes the LJ 12-6 potential with the well-depth @p epsilon_0
  !! and contact distance @p sigma_0.
  function lj_init_parameters(epsilon_0, sigma_0) result(this)
    real(dp), intent(in) :: epsilon_0, sigma_0
    type(lj_potential) :: this
    this%epsilon_0 = epsilon_0
    this%sigma_0 = sigma_0
  end function lj_init_parameters
  
  !> Initialize the module using parameters given by @p reader.
  function lj_init_wt_reader(reader) result(this)
    type(parameterizer), intent(in) :: reader
    type(lj_potential) :: this
    call getparameter(reader, 'lj_sigma_0', this%sigma_0)
    call getparameter(reader, 'lj_epsilon_0', this%epsilon_0)
  end function lj_init_wt_reader
  
  !> Write the module parameters using the output unit and format defined
  !! by @p writer.
  subroutine lj_writeparameters(this, writer)
    class(lj_potential), intent(in) :: this
    type(parameter_writer), intent(inout) :: writer
    call writecomment(writer, 'Lennard-Jones potential parameters.')
    call writeparameter(writer, 'lj_sigma_0', this%sigma_0)
    call writeparameter(writer, 'lj_epsilon_0', this%epsilon_0)
  end subroutine lj_writeparameters
  
  !> Returns the LJ 12-6 potential with the parameterization set with
  !! lj_init and a given internuclear vector @p rij.
  pure subroutine lj1(this, r, energy, overlap)
    class(lj_potential), intent(in) :: this
    real(dp), intent(in) :: r
    real(dp), intent(out) :: energy
    logical, intent(out), optional :: overlap
    energy = lj3(r, this%epsilon_0, this%sigma_0)
    if (present(overlap)) overlap = .false.
  end subroutine lj1
  
  !> Returns the LJ 12-6 potential value with given well-depth @p epsilon, 
  !! range parameter @p sigma and the internuclear vector @p rij from
  !! particle i to particle j. 
  pure function lj3(rij, epsilon, sigma)
    real(dp) :: lj3
    real(dp), intent(in) :: rij
    real(dp), intent(in) :: epsilon
    real(dp), intent(in) :: sigma
    lj3 = (sigma / rij)**6
    lj3 = lj3 * (lj3 - 1.)
    lj3 = 4. * epsilon * lj3
  end function lj3
  
  !> Returns the force exerted on the particle j by particle i when the
  !! interparticle vector from i to j is @p rij.
  pure function lj_force(this, rij)
    class(lj_potential), intent(in) :: this
    real(dp), intent(in) :: rij(3)
    real(dp) :: lj_force(3)
    real(dp) :: urij(3)
    real(dp) :: rijabs
    rijabs = sqrt(dot_product(rij, rij))
    urij = rij/rijabs
    lj_force = 24. * this%epsilon_0 * this%sigma_0**6 / rijabs**7 * &
         (2. * this%sigma_0**6 / rijabs**6 - 1.) * urij
  end function lj_force
  
end module m_lj
