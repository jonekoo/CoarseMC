module m_gb_interaction
  use m_gayberne, only: gayberne, gb_from_json
  use particle, only: pair_interaction, particledat
  use num_kind, only: dp
  use json_module
  use class_parameter_writer, only: parameter_writer, writeparameter
  use m_json_wrapper, only: get_parameter
  implicit none

  type, extends(pair_interaction) :: gb_interaction
     type(gayberne) :: pef
     real(dp) :: cutoff = 5.5
   contains
     procedure :: pair_potential => gb_pair_potential
     procedure :: pair_force => gb_pair_force
     procedure :: get_cutoff => gb_get_cutoff
     procedure :: to_json => gb_to_json
     procedure :: writeparameters => gb_writeparameters
  end type gb_interaction

  interface gb_interaction
     module procedure gb_interaction_from_json
  end interface gb_interaction
  
contains

  function gb_interaction_from_json(json_val) result(this)
    type(json_value), pointer, intent(in) :: json_val
    type(gb_interaction) :: this
    call get_parameter(json_val, 'r_cutoff', this%cutoff, error_lb=0._dp)
    this%pef = gayberne(json_val)
  end function gb_interaction_from_json

  subroutine gb_to_json(this, json_val)
    class(gb_interaction), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    type(json_value), pointer :: pef_json
    call json_add(json_val, 'type', 'gb_interaction')
    call json_add(json_val, 'r_cutoff', this%cutoff)
    call json_create_object(pef_json, 'potential')
    call this%pef%to_json(pef_json)
    call json_add(json_val, pef_json)
  end subroutine gb_to_json
  
  pure subroutine gb_pair_potential(this, particlei, particlej, rij, energy, err)
    class(gb_interaction), intent(in) :: this
    type(particledat), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    call this%pef%potential(particlei%orientation(), particlej%orientation(), &
         rij, energy, overlap)
    if (overlap) err = 1
  end subroutine gb_pair_potential

  pure function gb_get_cutoff(this) result(cutoff)
    class(gb_interaction), intent(in) :: this
    real(dp) :: cutoff
    cutoff = this%cutoff
  end function gb_get_cutoff

  pure function gb_pair_force(this, particlei, particlej, rij) result(f)
    class(gb_interaction), intent(in) :: this
    type(particledat), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp) :: f(3)
    f = this%pef%force(particlei%orientation(), particlej%orientation(), rij)
  end function gb_pair_force

  subroutine gb_writeparameters(this, writer)
    class(gb_interaction), intent(in) :: this
    type(parameter_writer), intent(inout) :: writer
    call writeparameter(writer, 'gb_r_cutoff', this%cutoff)
    call this%pef%writeparameters(writer)
  end subroutine gb_writeparameters
  
end module m_gb_interaction
