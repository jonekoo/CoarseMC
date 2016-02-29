!> Routines related to the Lennard-Jones interaction wrapper
!! lj_interaction.
module m_lj_interaction
  use m_lj, only: lj_potential, lj_from_json
  use m_particle, only: pair_interaction, particle
  use num_kind, only: dp
  use json_module
  use m_json_wrapper, only: get_parameter
  use m_rod, only: rod
  use m_point, only: point
  implicit none

  !> Wraps the Lennard-Jones potential as a pair_interaction.
  type, extends(pair_interaction) :: lj_interaction

     !> The Lennard-Jones potential.
     type(lj_potential) :: pef

     !> Cutoff radius of the interaction.
     real(dp) :: cutoff = 5.5
     
   contains

     !> Computes the energy of the interaction between two LJ particles.
     procedure :: pair_potential => lj_pair_potential

     !> Computes the force acting between two LJ particles.
     procedure :: pair_force => lj_pair_force

     !> Returns the cutoff radius of the interaction.
     procedure :: get_cutoff => lj_get_cutoff

     !> Serializes the interaction to JSON.
     procedure :: to_json => lj_to_json

     !> Produces a sample of the values of this interaction.
     procedure :: sample => lj_sample
  end type lj_interaction

  !> Contstuctor interface for the LJ interaction.
  interface lj_interaction
     module procedure lj_interaction_from_json
  end interface lj_interaction
  
contains

  !> Deserializes the lj_interaction from JSON value @p json_val.
  function lj_interaction_from_json(json_val) result(this)
    type(json_value), pointer, intent(in) :: json_val
    type(lj_interaction) :: this
    call get_parameter(json_val, 'r_cutoff', this%cutoff, error_lb=0._dp)
    this%pef = lj_potential(json_val)
  end function lj_interaction_from_json


  !> Serializes @p this interaction to JSON value @p json_val.
  subroutine lj_to_json(this, json_val)
    class(lj_interaction), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    call json_add(json_val, 'r_cutoff', this%cutoff)
    call this%pef%to_json(json_val)
  end subroutine lj_to_json


  !> Computes the @p energy of @p this interaction for the particle pair
  !! @p particlei, @p particlej. @p rij is the vector from @p particlei
  !! to @p particlej. If @p err = 1 the particles are considered to
  !! overlap.
  pure subroutine lj_pair_potential(this, particlei, particlej, rij, &
       energy, err)
    class(lj_interaction), intent(in) :: this
    class(particle), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    err = 0
    call this%pef%potential(norm2(rij), energy, overlap)
    if (overlap) err = 1
  end subroutine lj_pair_potential


  !> Returns the cutoff radius of @p this interaction.
  pure function lj_get_cutoff(this) result(cutoff)
    class(lj_interaction), intent(in) :: this
    real(dp) :: cutoff
    cutoff = this%cutoff
  end function lj_get_cutoff


  !> Returns the force acting on @p particlei by the @p particlej due to
  !! @p this interaction. @p rij is the vector from @p particlei to
  !! @p particlej.
  pure function lj_pair_force(this, particlei, particlej, rij) result(f)
    class(lj_interaction), intent(in) :: this
    class(particle), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp) :: f(3)
    f = this%pef%force(rij)
  end function lj_pair_force


  !> Produces a sample of the possible potential energies of @p this
  !! interaction with particle-to-particle distances r. The sample is
  !! written to JSON value @p json_val.
  subroutine lj_sample(this, json_val, r)
    class(lj_interaction), intent(in) :: this
    type(json_value), pointer :: json_val
    type(point) :: lj
    real(dp), intent(in) :: r(:)
    type(json_value), pointer :: child
    real(dp) :: energy
    integer :: err, i, j, k
    type(json_value), pointer :: r_json, energy_json

    call json_create_object(json_val, 'lj_interaction')
    call json_create_object(child, 'LJ-LJ')
    call json_create_array(r_json, 'x')
    call json_create_array(energy_json, 'energy')
    do k = 1, size(r)
       call this%pair_potential(lj, lj, [r(k), 0._dp, 0._dp], energy, err)
       if (err == 0) then
          call json_add(r_json, '', r(k))
          call json_add(energy_json, '', energy)
       end if
    end do
    call json_add(child, r_json)
    call json_add(child, energy_json)
    call json_add(json_val, child)
  end subroutine lj_sample
  

end module m_lj_interaction
