!> Routines related to the Gay-Berne interaction wrapper gb_interaction.
module m_gb_interaction
  use m_gayberne, only: gayberne, gb_from_json
  use m_point, only: pair_interaction, point
  use num_kind, only: dp
  use json_module
  use m_json_wrapper, only: get_parameter
  use m_rod, only: rod
  implicit none

  !> Wraps the Gay-Berne potential as a pair_interaction.
  type, extends(pair_interaction) :: gb_interaction

     !> The Gay-Berne potential
     type(gayberne) :: pef

     !> Cutoff distance of the interaction
     real(dp) :: cutoff = 5.5
     
   contains
     
     !> Computes the potential energy of the interaction of two rods.
     procedure :: pair_potential => gb_pair_potential

     !> Computes the force acting between two rods due to this
     !! interaction. 
     procedure :: pair_force => gb_pair_force

     !> Returns the cutoff of this interaction.
     procedure :: get_cutoff => gb_get_cutoff

     !> Serializes this interaction to JSON.
     procedure :: to_json => gb_to_json

     !> Produces a sample of the energies of the Gay-Berne interaction. 
     procedure :: sample => gb_sample
  end type gb_interaction

  interface gb_interaction
     module procedure gb_interaction_from_json
  end interface gb_interaction
  
contains

  !> Deserializes the interaction from the JSON value @p json_val.
  function gb_interaction_from_json(json_val) result(this)
    type(json_value), pointer, intent(in) :: json_val
    type(gb_interaction) :: this
    call get_parameter(json_val, 'r_cutoff', this%cutoff, error_lb=0._dp)
    this%pef = gayberne(json_val)
  end function gb_interaction_from_json

  !> Serializes @p this interaction to the JSON value @p json.
  subroutine gb_to_json(this, json_val)
    class(gb_interaction), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    call json_add(json_val, 'r_cutoff', this%cutoff)
    call this%pef%to_json(json_val)
  end subroutine gb_to_json

  
  !> Computes the potential @p energy of the Gay-Berne interaction
  !! acting between @p particlei and @p particlej. @p rij is the vector
  !! from @p particlei to @p particlej. If the particles overlap,
  !! @p err = 1. If particlej is not an instance of class rod,
  !! @p err = 78. If particlej is not an instance of class rod
  !! @p err = 77.
  !! 
  pure subroutine gb_pair_potential(this, particlei, particlej, rij, &
       energy, err)
    class(gb_interaction), intent(in) :: this
    class(point), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    err = 0
    select type (particlei)
    class is (rod)
       select type (particlej)
       class is (rod)
          call this%pef%potential(particlei%orientation(), &
               particlej%orientation(), rij, energy, overlap)
          if (overlap) err = 1
       class default
         err = 77
       end select
    class default
      err = 78
    end select
  end subroutine gb_pair_potential


  !> Returns the cutoff radius of @p this interaction.
  pure function gb_get_cutoff(this) result(cutoff)
    class(gb_interaction), intent(in) :: this
    real(dp) :: cutoff
    cutoff = this%cutoff
  end function gb_get_cutoff


  !> Returns the force acting on @p particlei by @p particlej due to
  !! @p this interaction between the particles. @p rij is the vector from
  !! particle i to particle j.
  pure function gb_pair_force(this, particlei, particlej, rij) result(f)
    class(gb_interaction), intent(in) :: this
    class(point), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp) :: f(3)
    select type (particlei)
    class is (rod)
       select type (particlej)
       class is (rod)
          f = this%pef%force(particlei%orientation(), &
               particlej%orientation(), rij)
       end select
    end select
  end function gb_pair_force


  !> Produces a sample of the possible potential energies with given
  !! distances between particles @p r. The sample is written to the JSON
  !! value @p json_val. The sample includes combinations of particles
  !! oriented in x, y, and z directions.
  subroutine gb_sample(this, json_val, r)
    class(gb_interaction), intent(in) :: this
    type(json_value), pointer :: json_val
    type(rod) :: rodsamples(3)
    real(dp), intent(in) :: r(:)
    type(json_value), pointer :: child
    real(dp) :: energy
    integer :: err, i, j, k
    type(json_value), pointer :: r_json, energy_json
    character(len=3, kind=CK), parameter :: descriptions(3) = &
         ['GBx', 'GBy', 'GBz']

    call rodsamples(1)%set_orientation([1._dp, 0._dp, 0._dp])
    call rodsamples(2)%set_orientation([0._dp, 1._dp, 0._dp])
    call rodsamples(3)%set_orientation([0._dp, 0._dp, 1._dp])
    call json_create_object(json_val, 'gb_interaction')
    do i = 1, size(rodsamples), 2
       do j = 1, i
          call json_create_object(child, descriptions(i) // '-' // &
               descriptions(j))
          call json_create_array(r_json, 'x')
          call json_create_array(energy_json, 'energy')
          do k = 1, size(r)
             call this%pair_potential(rodsamples(i), rodsamples(j), &
                  [r(k), 0._dp, 0._dp], energy, err)
             if (err == 0) then
                call json_add(r_json, '', r(k))
                call json_add(energy_json, '', energy)
             end if
          end do
          call json_add(child, r_json)
          call json_add(child, energy_json)
          call json_add(json_val, child)
       end do
    end do
  end subroutine gb_sample
  
end module m_gb_interaction
