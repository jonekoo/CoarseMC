!> Module responsible for defining the computation of the (uniaxial)
!! Gay-Berne potential for the interaction of two ellipsoidal molecules.
!!
!! @see J. G. Gay and B. J. Berne. J. Chem. Phys., 74(6):3316, 1981.
!! @see M. A. Bates and G. R. Luckhurst. J. Chem. Phys.,
!!      110(14):7087, 1999.
!!
module m_gayberne_interaction
  use nrtype, only: dp
  use class_parameterizer
  use class_parameter_writer
  use m_rod
  use m_stopping_internal_interaction
  use m_gayberne, only: gayberne_potential
  use json_module
  implicit none

  type, extends(stopping_internal_interaction) :: gayberne_interaction
     type(gayberne_potential) :: gb_potential
     !> Defines the GB distance when the particles overlap in
     !! gb_potential.
     real(dp) :: hardcore = 0.6_dp

     real(dp) :: cutoff = 5.5
     
   contains
     procedure :: rod_potential => gayberne_interaction_potential
     procedure :: writeparameters => gb_writeparameters
     procedure :: get_cutoff => gb_get_cutoff
     procedure :: from_json => gb_from_json
     procedure :: to_json => gb_to_json
  end type gayberne_interaction

  !> Initializes the type
  interface gayberne_interaction
     module procedure initparameterizer, construct_gb_interaction
  end interface gayberne_interaction
  
  !> The interface for the well-depth function in the GB potential.
  interface
     pure function gb_eps(this, ui, uj, urij)
       import gayberne_interaction, dp
       implicit none
       class(gayberne_interaction), intent(in) :: this 
       real(dp), intent(in) :: ui(3), uj(3), urij(3)
       real(dp) :: gb_eps
     end function gb_eps
  end interface

contains

  function construct_gb_interaction(gb_potential, hardcore, cutoff) result(this)
    type(gayberne_potential), intent(in) :: gb_potential
    real(dp), intent(in) :: hardcore, cutoff
    type(gayberne_interaction) :: this
    this%gb_potential = gb_potential
    this%hardcore = hardcore
    this%cutoff = cutoff
  end function construct_gb_interaction


  subroutine gb_from_json(this, json_val)
    class(gayberne_interaction), intent(inout) :: this
    type(json_value), pointer, intent(in) :: json_val
  end subroutine gb_from_json

  subroutine gb_to_json(this, json_val)
    class(gayberne_interaction), intent(in) :: this
    type(json_value), pointer :: json_val
    !! Write name of interaction
    call json_add(json_val, 'name', 'gayberne_interaction')
    !! write GB potential parameters
    !! Write extra parameters: hardcore, cutoff,...
    !! Write group name
  end subroutine gb_to_json

  
  !> Initializes the module using a parameterizer object.
  !! 
  !! @param[in] reader the parameterizer object which is responsible
  !! for getting the parameters for this module from some source, e.g.
  !! file. 
  !!
  function initparameterizer(reader) result(this)
    type(gayberne_interaction) :: this 
    type(parameterizer), intent(in) :: reader
    real(dp) :: gb_kappasigma
    real(dp) :: gb_kappaepsilon
    real(dp) :: gb_mu
    real(dp) :: gb_nu
    real(dp) :: gb_sigma_0
    real(dp) :: gb_epsilon_0
    call getparameter(reader, 'gb_kappa_sigma', gb_kappasigma)
    call getparameter(reader, 'gb_kappa_epsilon', gb_kappaepsilon)
    call getparameter(reader, 'gb_mu', gb_mu)
    call getparameter(reader, 'gb_nu', gb_nu)
    call getparameter(reader, 'gb_sigma_0', gb_sigma_0)
    call getparameter(reader, 'gb_epsilon_0', gb_epsilon_0)
    call getparameter(reader, 'gb_hardcore', this%hardcore)
    call getparameter(reader, 'gb_cutoff', this%cutoff)
    this%gb_potential = gayberne_potential(gb_kappasigma, gb_kappaepsilon, &
         gb_mu, gb_nu, gb_sigma_0, gb_epsilon_0) 
  end function initparameterizer

  !> Write the parameters of this module to the output unit and format
  !! defined by @p writer.
  subroutine gb_writeparameters(this, writer)
    class(gayberne_interaction), intent(in) :: this
    type(parameter_writer), intent(inout) :: writer
    call writecomment(writer, 'Gay-Berne potential parameters')
    call writeparameter(writer, 'gb_kappa_sigma', this%gb_potential%kappasigma)
    call writeparameter(writer, 'gb_kappa_epsilon', &
         this%gb_potential%kappaepsilon)
    call writeparameter(writer, 'gb_mu', this%gb_potential%mu)
    call writeparameter(writer, 'gb_nu', this%gb_potential%nu)
    call writeparameter(writer, 'gb_sigma_0', this%gb_potential%sigma0)
    call writeparameter(writer, 'gb_epsilon_0', this%gb_potential%epsilon0)    
    call writeparameter(writer, 'gb_hardcore', this%hardcore)
    call writeparameter(writer, 'gb_cutoff', this%cutoff)
  end subroutine gb_writeparameters

  subroutine gayberne_interaction_potential(this, rod_i, rod_j, &
       res, err)
    class(gayberne_interaction), intent(in) :: this
    type(rod), intent(in) :: rod_i(:), rod_j(:)
    real(dp), intent(out) :: res
    integer, intent(out) :: err
    logical :: overlap
    real(dp) :: s
    integer :: i, j
    err = 0
    res = 0.
    !! Could optimize this by precomputing
    !! 1) pairs inside cutoff distance
    !! 2) overlapping particles
    !! As well as making direct references to orientation and position.
    do i = 1, size(rod_i)
       do j = 1, size(rod_j)
          associate(r => rod_j(j)%position() - rod_i(i)%position())
            if (norm2(r) < this%cutoff) then 
               call this%gb_potential%potential(rod_i(i)%orientation(), &
                    rod_j(j)%orientation(), r, s, overlap)
               if (overlap) err = -1
               res = res + s   
            end if
          end associate
       end do
    end do
  end subroutine gayberne_interaction_potential

  
  function gb_get_cutoff(this)
    class(gayberne_interaction), intent(in) :: this
    real(dp) :: gb_get_cutoff
    gb_get_cutoff = this%cutoff
  end function gb_get_cutoff
  
  
end module m_gayberne_interaction


