module m_gblj_interaction
  use m_gblj
  use m_particle, only: pair_interaction, particle
  use num_kind, only: dp
  use json_module
  use class_parameter_writer, only: parameter_writer, writeparameter
  use m_json_wrapper, only: get_parameter
  use m_rod, only: rod
  use m_point, only: point
  implicit none

  type, extends(pair_interaction) :: gblj_interaction
     type(gblj_potential) :: pef
     real(dp) :: cutoff = 5.5
   contains
     procedure :: pair_potential => gblj_pair_potential
     procedure :: pair_force => gblj_pair_force
     procedure :: get_cutoff => gblj_get_cutoff
     procedure :: to_json => gblj_to_json
  end type gblj_interaction

  interface gblj_interaction
     module procedure gblj_interaction_from_json
  end interface gblj_interaction
  
contains

  function gblj_interaction_from_json(json_val) result(this)
    type(json_value), pointer, intent(in) :: json_val
    type(gblj_interaction) :: this
    call get_parameter(json_val, 'r_cutoff', this%cutoff, error_lb=0._dp)
    this%pef = gblj_potential(json_val)
  end function gblj_interaction_from_json

  subroutine gblj_to_json(this, json_val)
    class(gblj_interaction), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    call json_add(json_val, 'r_cutoff', this%cutoff)
    call this%pef%to_json(json_val)
  end subroutine gblj_to_json
  
  pure subroutine gblj_pair_potential(this, particlei, particlej, rij, &
       energy, err)
    class(gblj_interaction), intent(in) :: this
    class(particle), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    err = 0
    select type (particlei)
    class is (rod)
          call this%pef%potential(particlei%orientation(), rij, energy, overlap)
          if (overlap) err = 1
    class default
       select type (particlej)
       class is (rod)
          call this%pef%potential(particlej%orientation(), -rij, energy, &
               overlap)
          if (overlap) err = 1
       class default
         err = 79
       end select
    end select
  end subroutine gblj_pair_potential

  pure function gblj_get_cutoff(this) result(cutoff)
    class(gblj_interaction), intent(in) :: this
    real(dp) :: cutoff
    cutoff = this%cutoff
  end function gblj_get_cutoff

  pure function gblj_pair_force(this, particlei, particlej, rij) result(f)
    class(gblj_interaction), intent(in) :: this
    class(particle), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp) :: f(3)
    f = 0._dp
    select type (particlei)
    type is (rod)
       select type (particlej)
       type is (point)
          f = this%pef%force(particlei%orientation(), rij)
       end select
    type is (point)
       select type (particlej)
       type is (rod)
          f = this%pef%force(particlej%orientation(), -rij)
       end select
    end select
  end function gblj_pair_force

end module m_gblj_interaction
