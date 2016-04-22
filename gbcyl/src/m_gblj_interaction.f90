!> Gay-Berne - Lennard-Jones interaction gblj_interaction between a
!! rod-like and a point-like particle and related procedures.
module m_gblj_interaction
  use m_gblj
  use m_particle, only: pair_interaction, point
  use num_kind, only: dp
  use json_module
  use m_json_wrapper, only: get_parameter
  use m_rod, only: rod
  implicit none

  !> Wraps the GB-LJ potential as a pair_interaction.
  type, extends(pair_interaction) :: gblj_interaction

     !> The GB-LJ potential.
     type(gblj_potential) :: pef

     !> The cutoff radius of the interaction.
     real(dp) :: cutoff = 5.5

   contains

     !> Computes the GB-LJ interaction energy.
     procedure :: pair_potential => gblj_pair_potential

     !> Computes the force acting between a pair of particles.
     procedure :: pair_force => gblj_pair_force

     !> Returns the cutoff radius of this interaction.
     procedure :: get_cutoff => gblj_get_cutoff

     !> Serializes this interaction to JSON.
     procedure :: to_json => gblj_to_json

     !> Returns a sample of the possible potential energies of the
     !! interaction.
     procedure :: sample => gblj_sample
  end type gblj_interaction

  !> Constructor interface.
  interface gblj_interaction
     module procedure gblj_interaction_from_json
  end interface gblj_interaction
  
contains

  !> Deserializes the GB-LJ interaction from the JSON value @p json_val. 
  function gblj_interaction_from_json(json_val) result(this)
    type(json_value), pointer, intent(in) :: json_val
    type(gblj_interaction) :: this
    call get_parameter(json_val, 'r_cutoff', this%cutoff, error_lb=0._dp)
    this%pef = gblj_potential(json_val)
  end function gblj_interaction_from_json


  !> Serializes @p this interactions to JSON value @p json_val.
  subroutine gblj_to_json(this, json_val)
    class(gblj_interaction), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    call json_add(json_val, 'r_cutoff', this%cutoff)
    call this%pef%to_json(json_val)
  end subroutine gblj_to_json
  
  !> Computes the @p energy of @p this interaction between @p particlei
  !! and @p particlej. @p rij is the vector from @p particlei to 
  !! @p particlej. @p err = 1 if the particles overlap. @p err = 79 if
  !! neither @p particlei or @p particlej is a rod.
  pure subroutine gblj_pair_potential(this, particlei, particlej, rij, &
       energy, err)
    class(gblj_interaction), intent(in) :: this
    class(point), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    err = 0
    select type (particlei)
    class is (rod)
          call this%pef%potential(particlei%orientation(), rij, energy, &
               overlap)
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


  !> Returns the cutoff radius of @p this interaction.
  pure function gblj_get_cutoff(this) result(cutoff)
    class(gblj_interaction), intent(in) :: this
    real(dp) :: cutoff
    cutoff = this%cutoff
  end function gblj_get_cutoff


  !> Returns the force acting on @p particlei by @p particlej due to
  !! @p this interaction between them. @p rij is the vector from
  !! @p particlei to @p particlej.
  pure function gblj_pair_force(this, particlei, particlej, rij) result(f)
    class(gblj_interaction), intent(in) :: this
    class(point), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp) :: f(3)
    f = 0._dp
    select type (particlei)
    class is (rod)
       f = this%pef%force(particlei%orientation(), rij)
    class default
       select type (particlej)
       class is (rod)
          f = this%pef%force(particlej%orientation(), -rij)
       end select
    end select
  end function gblj_pair_force


  !> Produces a sample of the possible potential energies with @p this
  !! interaction with rod-to-point distances @p r. The results are
  !! stored to JSON value @p json_val. Datasets are created for rods
  !! oriented along x and y axes.
  subroutine gblj_sample(this, json_val, r)
    class(gblj_interaction), intent(in) :: this
    type(json_value), pointer :: json_val
    type(rod) :: rodsamples(2)
    type(point) :: lj
    real(dp), intent(in) :: r(:)
    type(json_value), pointer :: child
    real(dp) :: energy
    integer :: err, i, k
    type(json_value), pointer :: r_json, energy_json
    character(len=3, kind=CK), parameter :: descriptions(2) = &
         ['GBx', 'GBy']

    call rodsamples(1)%set_orientation([1._dp, 0._dp, 0._dp])
    call rodsamples(2)%set_orientation([0._dp, 1._dp, 0._dp])
    call json_create_object(json_val, 'gblj_interaction')
    do i = 1, size(rodsamples)   
       call json_create_object(child, descriptions(i) // '-LJ')
       call json_create_array(r_json, 'x')
       call json_create_array(energy_json, 'energy')
       do k = 1, size(r)
          call this%pair_potential(rodsamples(i), lj, &
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
  end subroutine gblj_sample

  
end module m_gblj_interaction
