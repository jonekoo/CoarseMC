!> Routines related to the Gay-Berne interaction wrapper gb_interaction.
module m_dgb_interaction
  use m_gayberne, only: gayberne, gb_from_json
  use m_point, only: pair_interaction, point
  use num_kind, only: dp
  use json_module
  use m_json_wrapper, only: get_parameter
  use m_rod, only: rod
  use m_doublerod, only: doublerod
  implicit none

  !> Wraps the Gay-Berne potential as a pair_interaction.
  type, extends(pair_interaction) :: dgb_interaction

     !> The Gay-Berne potential
     type(gayberne) :: pef

     !> Cutoff distance of the interaction
     real(dp) :: cutoff = 5.5 !+length of gb 
     
   contains
     
     !> Computes the potential energy of the interaction of two rods.
     procedure :: pair_potential => dgb_pair_potential

     !> Computes the force acting between two rods due to this
     !! interaction. 
     procedure :: pair_force => dgb_pair_force

     !> Returns the cutoff of this interaction.
     procedure :: get_cutoff => dgb_get_cutoff

     !> Serializes this interaction to JSON.
     procedure :: to_json => dgb_to_json

     !> Produces a sample of the energies of the Gay-Berne interaction. 
     procedure :: sample => dgb_sample
  end type dgb_interaction

  interface dgb_interaction
     module procedure dgb_interaction_from_json
  end interface dgb_interaction
  
contains

  !> Deserializes the interaction from the JSON value @p json_val.
  function dgb_interaction_from_json(json_val) result(this)
    type(json_value), pointer, intent(in) :: json_val
    type(dgb_interaction) :: this
    call get_parameter(json_val, 'r_cutoff', this%cutoff, error_lb=0._dp)
    this%pef = gayberne(json_val)

    !lisää tarvittavat parametrit tähän!
  end function dgb_interaction_from_json

  !> Serializes @p this interaction to the JSON value @p json.
  subroutine dgb_to_json(this, json_val)
    class(dgb_interaction), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    call json_add(json_val, 'r_cutoff', this%cutoff)
    call this%pef%to_json(json_val)
  end subroutine dgb_to_json

  
  !> Computes the potential @p energy of the Gay-Berne interaction
  !! acting between @p particlei and @p particlej. @p rij is the vector
  !! from @p particlei to @p particlej. If the particles overlap,
  !! @p err = 1. If particlej is not an instance of class rod,
  !! @p err = 78. If particlej is not an instance of class rod
  !! @p err = 77.
  !! 
  pure subroutine dgb_pair_potential(this, particlei, particlej, rij, &
       energy, err)
    class(dgb_interaction), intent(in) :: this
    class(point), intent(in) :: particlei, particlej
    real(dp) :: gblength, potE1, potE2, potE3, potE4 
    integer :: i
    real(dp), intent(in) :: rij(3)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    real(dp), dimension(3) :: ui, uj, vi, vj, dist1, dist2, dist3, dist4
    err = 0
    select type (particlei)
    class is (doublerod)
      select type (particlej)
      class is (doublerod)
        gblength = this%pef%kappasigma !! from gayberne, I presume

        !! this should be made a bit more optimised: now local copies
        !! of orientation vectors are created.
        ui = [particlei%ux, particlei%uy, particlei%uz]
        vi = [particlei%vx, particlei%vy, particlei%vz]
        uj = [particlej%ux, particlej%uy, particlej%uz]
        vj = [particlej%vx, particlej%vy, particlej%vz]

        do i = 1, 3
          dist1(i) = rij(i) + gblength/2*(uj(i) - vi(i))
          dist2(i) = rij(i) + gblength/2*(vj(i) - vi(i))
          dist3(i) = rij(i) + gblength/2*(vj(i) - ui(i))
          dist4(i) = rij(i) + gblength/2*(uj(i) - ui(i))
        end do

        call this%pef%potential(vi, uj, dist1, potE1, overlap)
        call this%pef%potential(vi, vj, dist2, potE2, overlap)
        call this%pef%potential(ui, vj, dist3, potE3, overlap)
        call this%pef%potential(ui, uj, dist4, potE4, overlap)
        energy = potE1 + potE2 + potE3 + potE4
        if (overlap) err = 1

      class default
        err = 177
      end select
    class default
      err = 178
    end select
  end subroutine dgb_pair_potential


  !> Returns the cutoff radius of @p this interaction.
  pure function dgb_get_cutoff(this) result(cutoff)
    class(dgb_interaction), intent(in) :: this
    real(dp) :: cutoff
    cutoff = this%cutoff
  end function dgb_get_cutoff


  !> Returns the force acting on @p particlei by @p particlej due to
  !! @p this interaction between the particles. @p rij is the vector from
  !! particle i to particle j.
  !! OBS! HUOM! NB! THIS ISN'T CONVERTED YET FOR DGB MOLECULES!

  pure function dgb_pair_force(this, particlei, particlej, rij) result(f)
    class(dgb_interaction), intent(in) :: this
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
  end function dgb_pair_force


  !> Produces a sample of the possible potential energies with given
  !! distances between particles @p r. The sample is written to the JSON
  !! value @p json_val. The sample includes combinations of particles
  !! oriented in x, y, and z directions.
  subroutine dgb_sample(this, json_val, r)
    class(dgb_interaction), intent(in) :: this
    type(json_value), pointer :: json_val
    type(doublerod) :: doublerodsamples(3)
    real(dp), intent(in) :: r(:)
    type(json_value), pointer :: child
    real(dp) :: energy
    integer :: err, i, j, k
    type(json_value), pointer :: r_json, energy_json
    character(len=3, kind=CK), parameter :: descriptions(3) = &
         ['DGBx', 'DGBy', 'DGBz']

    call doublerodsamples(1)%set_orientation_u([1._dp, 0._dp, 0._dp])
    call doublerodsamples(1)%set_orientation_v([-1._dp, 0._dp, 0._dp])
    call doublerodsamples(2)%set_orientation_u([0._dp, 1._dp, 0._dp])
    call doublerodsamples(2)%set_orientation_v([0._dp, -1._dp, 0._dp])
    call doublerodsamples(3)%set_orientation_u([0._dp, 0._dp, 1._dp])
    call doublerodsamples(3)%set_orientation_v([0._dp, 0._dp, -1._dp])
    call json_create_object(json_val, 'dgb_interaction')
    do i = 1, size(doublerodsamples), 2
       do j = 1, i
          call json_create_object(child, descriptions(i) // '-' // &
               descriptions(j))
          call json_create_array(r_json, 'x')
          call json_create_array(energy_json, 'energy')
          do k = 1, size(r)
             call this%pair_potential(doublerodsamples(i), doublerodsamples(j), &
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
  end subroutine dgb_sample
  
end module m_dgb_interaction
