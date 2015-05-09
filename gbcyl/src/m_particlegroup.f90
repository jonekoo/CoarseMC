module m_particlegroup
  use nrtype, only: dp
  use m_particle
  use class_poly_box
  use m_external_interaction, only: external_interaction
  use m_internal_interaction, only: internal_interaction
  include 'rng.inc'
  implicit none


  type, abstract :: energy_rule
   contains
     procedure(accept_by_de), deferred :: accept
  end type energy_rule

  interface
     pure subroutine accept_by_de(this, de, genstate, accepted)
       import energy_rule, dp, rngstate
       class(energy_rule), intent(inout) :: this
       real(dp), intent(in) :: de
       type(rngstate), intent(inout) :: genstate
       logical, intent(out) :: accepted
     end subroutine accept_by_de
  end interface
  
  type, abstract :: particlegroup
   contains
     procedure(move_with_rule), deferred  :: move
     procedure(modify_particles), deferred :: modify
     procedure(get_particle_count), deferred :: get_size
     procedure(get_particles), deferred :: get_particles 
     procedure(map_object), deferred :: map
     procedure(uncut_interact), private, deferred :: uncut_interact
     procedure(cut_interact), private, deferred :: cut_interact
     procedure(cut_pair_interact), private, deferred :: cut_pair_interact
     procedure(uncut_pair_interact), private, deferred :: uncut_pair_interact
     generic :: interact => uncut_interact, cut_interact
     generic :: pair_interact => cut_pair_interact, uncut_pair_interact
     procedure(add_internal_interaction), deferred :: add_internal_interaction
     procedure(add_external_interaction), deferred :: add_external_interaction
     generic :: add_interaction => add_internal_interaction, &
          add_external_interaction
  end type particlegroup

  interface
     subroutine add_internal_interaction(this, ia)
       import particlegroup, internal_interaction
       class(particlegroup), intent(inout) :: this
       class(internal_interaction), pointer, intent(in) :: ia
     end subroutine add_internal_interaction
  end interface
  
  interface
     subroutine add_external_interaction(this, ia)
       import particlegroup, external_interaction
       class(particlegroup), intent(inout) :: this
       class(external_interaction), pointer, intent(in) :: ia
     end subroutine add_external_interaction
  end interface
  
  interface
     subroutine get_particles(this, particles)
       import particlegroup, particle
       class(particlegroup), intent(in) :: this
       class(particle), allocatable, intent(out) :: particles(:)
     end subroutine get_particles
  end interface
  
  interface
     subroutine move_with_rule(this, rule, genstates, de)
       import particlegroup, energy_rule, rngstate, dp
       class(particlegroup), intent(inout) :: this
       class(energy_rule), intent(inout) :: rule
       type(rngstate), intent(inout) :: genstates(:)
       real(dp), intent(out) :: de
     end subroutine move_with_rule
  end interface
  
  interface
     subroutine particle_modifier(p)
       import particle
       class(particle), intent(inout) :: p
     end subroutine particle_modifier
  end interface

  interface
     subroutine modify_particles(this, proc)
       import particlegroup, particle_modifier
       class(particlegroup), intent(inout) :: this
       procedure(particle_modifier) :: proc
     end subroutine modify_particles
  end interface

  interface
     function get_particle_count(this)
       import particlegroup
       class(particlegroup), intent(in) :: this
       integer :: get_particle_count
     end function get_particle_count
  end interface


  type, abstract :: particle_handler
   contains
     procedure(handler_subr), deferred :: handler_subr
  end type particle_handler

  interface
     subroutine handler_subr(this, p)
       import particle_handler, particle
       class(particle_handler), intent(inout) :: this
       class(particle), intent(inout) :: p
     end subroutine handler_subr
  end interface

  interface
     subroutine map_object(this, obj)
       import particlegroup, particle_handler
       class(particlegroup), intent(inout) :: this
       class(particle_handler), intent(inout) :: obj
     end subroutine map_object
  end interface

  interface
     subroutine uncut_interact(this, ia, res, err)
       import particlegroup, external_interaction, dp
       class(particlegroup), intent(inout) :: this
       class(external_interaction), intent(inout) :: ia
       real(dp), intent(out) :: res
       integer, intent(out) :: err
     end subroutine uncut_interact
  end interface
  
  interface
     subroutine cut_interact(this, ia, position, cutoff, res, err)
       import particlegroup, external_interaction, dp
       class(particlegroup), intent(inout) :: this
       class(external_interaction), intent(inout) :: ia
       real(dp), intent(in) :: position(3)
       real(dp), intent(in) :: cutoff
       real(dp), intent(out) :: res
       integer, intent(out) :: err
     end subroutine cut_interact
  end interface

  interface
     subroutine cut_pair_interact(this, ia, cutoff, res, err)
       import particlegroup, internal_interaction, dp
       class(particlegroup), intent(inout) :: this
       class(internal_interaction), intent(inout) :: ia
       real(dp), intent(in) :: cutoff
       real(dp), intent(out) :: res
       integer, intent(out) :: err
     end subroutine cut_pair_interact
  end interface
  
  interface
     subroutine uncut_pair_interact(this, ia, res, err)
       import particlegroup, internal_interaction, dp
       class(particlegroup), intent(inout) :: this
       class(internal_interaction), intent(inout) :: ia
       real(dp), intent(out) :: res
       integer, intent(out) :: err
     end subroutine uncut_pair_interact
  end interface
  
  type, abstract :: particle_subroutine
   contains
     procedure(particle_subr), deferred :: subr
  end type particle_subroutine

  interface
     subroutine particle_subr(this, p)
       import particle_subroutine, particle
       class(particle_subroutine) :: this
       class(particle), intent(inout) :: p
     end subroutine particle_subr  
  end interface
    
  
end module m_particlegroup
