
  type, abstract :: MIXINPARAM_interaction
   contains
     procedure(MIXINPARAM_interaction_value), deferred :: value
  end type MIXINPARAM_interaction

  type, extends(MIXINPARAM) :: particlemixin
     class(MIXINPARAM_interaction), pointer :: interactions(:) => null()
   contains
     procedure :: energy => MIXINPARAM_energy
  end type particlemixin

  interface
     subroutine MIXINPARAM_interaction_value(this, the_MIXINPARAM, res, err)
       import MIXINPARAM_interaction, particlemixin, REAL64
       class(MIXINPARAM_interaction), intent(in) :: this
       type(particlemixin), intent(in) :: the_MIXINPARAM
       real(REAL64), intent(out) :: res
       integer, intent(out) :: err
     end subroutine MIXINPARAM_interaction_value
  end interface
  
