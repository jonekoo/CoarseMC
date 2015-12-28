module m_lj_interaction
  use m_lj, only: lj_potential, lj_from_json
  use particle, only: pair_interaction, particledat
  use num_kind, only: dp
  use json_module
  use m_json_wrapper, only: get_parameter
  use m_rod, only: rod
  use m_point, only: point
  implicit none

  type, extends(pair_interaction) :: lj_interaction
     type(lj_potential) :: pef
     real(dp) :: cutoff = 5.5
   contains
     procedure :: pair_potential => lj_pair_potential
     procedure :: pair_force => lj_pair_force
     procedure :: get_cutoff => lj_get_cutoff
     procedure :: to_json => lj_to_json
  end type lj_interaction

  interface lj_interaction
     module procedure lj_interaction_from_json
  end interface lj_interaction
  
contains

  function lj_interaction_from_json(json_val) result(this)
    type(json_value), pointer, intent(in) :: json_val
    type(lj_interaction) :: this
    call get_parameter(json_val, 'r_cutoff', this%cutoff, error_lb=0._dp)
    this%pef = lj_potential(json_val)
  end function lj_interaction_from_json

  subroutine lj_to_json(this, json_val)
    class(lj_interaction), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    call json_add(json_val, 'r_cutoff', this%cutoff)
    call this%pef%to_json(json_val)
  end subroutine lj_to_json
  
  pure subroutine lj_pair_potential(this, particlei, particlej, rij, &
       energy, err)
    class(lj_interaction), intent(in) :: this
    class(particledat), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    err = 0
    select type (particlei)
    type is (point)
       select type (particlej)
       type is (point)
          call this%pef%potential(norm2(rij), energy, overlap)
          if (overlap) err = 1
       class default
         err = 33
       end select
    class default
      err = 34
    end select
  end subroutine lj_pair_potential

  pure function lj_get_cutoff(this) result(cutoff)
    class(lj_interaction), intent(in) :: this
    real(dp) :: cutoff
    cutoff = this%cutoff
  end function lj_get_cutoff

  pure function lj_pair_force(this, particlei, particlej, rij) result(f)
    class(lj_interaction), intent(in) :: this
    class(particledat), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp) :: f(3)
    select type (particlei)
    type is (point)
       select type (particlej)
       type is (point)
          f = this%pef%force(rij)
       end select
    end select
  end function lj_pair_force

end module m_lj_interaction
