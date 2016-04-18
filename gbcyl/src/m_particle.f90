!> Contains the base type for all particles and related routines.
module m_particle
  use iso_fortran_env, only: error_unit, output_unit, dp => REAL64
  use utils, only: fmt_char_dp, acceptchange
  use class_poly_box, only: poly_box, minimage
  use mt_stream, only: rng=>genrand_double1_s, rngstate=>mt_state
  use particle_mover, only: transmove, rotate
  use json_module, only: json_value, json_add, CK, json_add, &
       json_create_array
  use m_json_wrapper, only: get_parameter
  implicit none

  !> The base type for all particles.
  type particle
     real(dp) :: x = 0._dp
     real(dp) :: y = 0._dp
     real(dp) :: z = 0._dp
   contains
     !> Returns the type of the particle as a string.
     procedure, nopass :: typestr => particle_typestr
     !> Returns a description of the coordinates of this particle.
     procedure, nopass :: description => particle_description
     !> Downcasts and assigns to particle. Overridden by subtypes.
     procedure :: downcast_assign => particle_downcast_assign
     !> Assignment operator implementation.
     procedure :: particle_assign
     generic :: assignment(=) => particle_assign
     
     !> Equality operator implementation. 
     procedure :: particle_equals
     generic :: operator(==) => particle_equals

     !> Deserializes the particle from string.
     procedure :: from_str => particle_from_str
     !> Deserializes the particle from JSON.
     procedure :: from_json => particle_from_json
     !> Serializes particle coordinates to JSON.
     procedure :: coordinates_to_json => particle_coordinates_to_json
     !> Serializes particle to the given Fortran output unit.
     procedure :: to_unit => particle_to_unit
     !> Returns the particle position.
     procedure :: position => position
     !> Sets the particle position.
     procedure :: set_position => particle_set_position
     !> Returns the potential energy of the with given neighbours and
     !! interactions.
     procedure :: energy => particle_potential
     procedure :: nvt_update => moveparticle_2
     !> Trial move generator. 
     procedure :: move => particle_move
  end type particle

  !> The base type for all pairwise-computed particle-particle
  !! interactions.
  type, abstract :: pair_interaction
   contains
     procedure(pair_potential), deferred :: pair_potential
     procedure(pair_force), deferred :: pair_force
     procedure(get_cutoff), deferred :: get_cutoff
     procedure(to_json), deferred :: to_json
  end type pair_interaction

  
  abstract interface

     !> Returns the pairwise-computed potential @p this for @p particlei
     !! and @p particlej.
     !!
     !! @param this the interaction.
     !! @param particlei, particlej the interacting particles.
     !! @param rij is the vector from particlei to particlej.
     !! @param energy the potential energy of the interaction.
     !! @param err non-zero if e.g. the particles overlap.
     !!
     pure subroutine pair_potential(this, particlei, particlej, rij, &
          energy, err)
       import
       class(pair_interaction), intent(in) :: this
       class(particle), intent(in) :: particlei, particlej
       real(dp), intent(in) :: rij(3)
       real(dp), intent(out) :: energy
       integer, intent(out) :: err
     end subroutine pair_potential

     
     !> Returns the force acting on @p particlei by the @p particlej
     !! resulting from interaction @p this.
     !!
     !! @param this the interaction.
     !! @param particlei, particlej the interacting particles.
     !! @param rij is the vector from particlei to particlej.
     !! @param energy the potential energy of the interaction.
     !! @param err non-zero if e.g. the particles overlap.
     !!
     pure function pair_force(this, particlei, particlej, rij) result(f)
       import
       class(pair_interaction), intent(in) :: this
       class(particle), intent(in) :: particlei, particlej
       real(dp), intent(in) :: rij(3)
       real(dp) :: f(3)
     end function pair_force


     !> Returns the cutoff radius of the @p this interaction. 
     pure function get_cutoff(this) result(res)
       import pair_interaction, dp
       class(pair_interaction), intent(in) :: this
       real(dp) :: res
     end function get_cutoff


     !> Serializes @p this interaction to JSON @p json_val.
     subroutine to_json(this, json_val)
       import pair_interaction, json_value
       class(pair_interaction), intent(in) :: this
       type(json_value), pointer, intent(inout) :: json_val
     end subroutine to_json
     
  end interface


  !> Wrapper for creating arrays of pair_interactions with different
  !! types.
  type pair_interaction_ptr
     class(pair_interaction), pointer :: ptr => null()
  end type pair_interaction_ptr


  !> Base type for all single-particle interactions.
  type, abstract :: single_interaction
   contains
     procedure(single_potential), deferred :: potential
     procedure(single_force), deferred :: force
     procedure(single_to_json), deferred :: to_json
  end type single_interaction
  
  abstract interface

     
     !> Returns the potential energy of @p particlei due to @p this
     !! interaction.
     !!
     !! @param this the interaction
     !! @param particlei the particle to which the energy is computed.
     !! @param simbox the simulation box.
     !! @param energy the potential energy of the particle.
     !! @param err non-zero if an error occurs.
     !! 
     pure subroutine single_potential(this, particlei, simbox, energy, err)
       import
       class(single_interaction), intent(in) :: this
       class(particle), intent(in) :: particlei
       type(poly_box), intent(in) :: simbox
       real(dp), intent(out) :: energy
       integer, intent(out) :: err
     end subroutine single_potential

     
     !> Returns the force acting on @p particlei due to @p this
     !! interaction.
     !!
     !! @param this interaction which causes the force.
     !! @param particlei the particle which the force acts on.
     !! @param simbox the simulation box in which the particle resides.
     !!
     pure function single_force(this, particlei, simbox) result(f)
       import
       class(single_interaction), intent(in) :: this
       class(particle), intent(in) :: particlei
       type(poly_box), intent(in) :: simbox
       real(dp) :: f(3)
     end function single_force


     !> Serializes @p this interaction to JSON value @p json_val.
     subroutine single_to_json(this, json_val)
       import single_interaction, json_value
       class(single_interaction), intent(in) :: this
       type(json_value), pointer, intent(inout) :: json_val
     end subroutine single_to_json
     
  end interface

  
  !> Wrapper to create arrays of single interactions of different types.
  type single_interaction_ptr
     class(single_interaction), pointer :: ptr => null()
  end type single_interaction_ptr

  
  !> Wrapper to be used for arrays of arrays of particles of different
  !! types.
  type particlearray_wrapper
     !> The particles.
     class(particle), allocatable :: arr(:)
     !> A mask which can be used to exclude particles.
     logical, allocatable :: mask(:)
     !> Number of particles to be considered in the array.
     integer :: n = 0
   contains
     procedure :: wrapper_assign
     generic :: assignment(=) => wrapper_assign
     procedure :: delete => wrapper_manual_delete
     final :: wrapper_delete
  end type particlearray_wrapper


contains

  !> Deserializes @p this particle from string @p str.
  !! If successful, ios==0.
  subroutine particle_from_str(this, str, ios)
    class(particle), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z
  end subroutine particle_from_str

  
  !> Returns the type of the particle as a string @p str.
  subroutine particle_typestr(str)
    character(len=:), allocatable, intent(out) :: str
    str = "particle"
  end subroutine particle_typestr


  !> Returns the description of particle coordinates as an array of
  !! strings. 
  subroutine particle_description(descr)
    character(kind=CK, len=3), allocatable, intent(inout) :: descr(:)
    descr = ["x  ", "y  ", "z  "]
  end subroutine particle_description

  
  !> Downcast and assign @p a_particle to @p this. If type of a_particle
  !! is different from @p this, @p err = 3.
  pure subroutine particle_downcast_assign(this, a_particle, err)
    class(particle), intent(inout) :: this
    class(particle), intent(in) :: a_particle
    integer, intent(out), optional :: err
    select type (a_particle)
    type is (particle)
       this = a_particle
       if (present(err)) err = 0
    class default
       if (present(err)) err = 3
    end select
  end subroutine particle_downcast_assign


  !> Assigns @p another to @p this.
  pure subroutine particle_assign(this, another)
    class(particle), intent(inout) :: this
    type(particle), intent(in) :: another
    this%x = another%x
    this%y = another%y
    this%z = another%z
  end subroutine particle_assign


  !> Returns .true. if @p another equals @p this.
  elemental function particle_equals(this, another) result(res)
    class(particle), intent(in) :: this
    type(particle), intent(in) :: another
    logical :: res
    res = (this%x == another%x) .and. (this%y == another%y) .and. &
         (this%z == another%z)
  end function particle_equals


  !> Translates @p this particle randomly. @p genstate is the random
  !! number generator state used to create the move.
  pure subroutine particle_move(this, genstate)    
    class(particle), intent(inout) :: this
    type(rngstate), intent(inout) :: genstate
    real(dp) :: xn, yn, zn
    call transmove(this%x, this%y, this%z, xn, yn, zn, genstate)
    this%x = xn
    this%y = yn
    this%z = zn
  end subroutine particle_move


  !> Serializes the coordinates of @p this particle to JSON value
  !! @p json_val.
  subroutine particle_coordinates_to_json(this, json_val)
    class(particle), intent(in) :: this
    type(json_value), pointer :: json_val
    call json_create_array(json_val, '')
    call json_add(json_val, '', this%x)
    call json_add(json_val, '', this%y)
    call json_add(json_val, '', this%z)
  end subroutine particle_coordinates_to_json


  !> Serializes the coordinates of @p this particle to output @p unit. 
  subroutine particle_to_unit(this, unit)
    class(particle), intent(in) :: this
    integer, intent(in) :: unit
    write(unit, *) 'particle', this%x, this%y, this%z
  end subroutine particle_to_unit


  !> Deserializes @p this particle from JSON value @p json_val.
  subroutine particle_from_json(this, json_val)
    class(particle), intent(inout) :: this
    type(json_value), pointer, intent(in) :: json_val
    call get_parameter(json_val, '[1]', this%x)
    call get_parameter(json_val, '[2]', this%y)
    call get_parameter(json_val, '[3]', this%z)
  end subroutine particle_from_json


  !> Serializes an array of @p particles to JSON value @p json_val.
  subroutine particlearray_to_json(json_val, particles)
    type(json_value), pointer, intent(inout) :: json_val
    class(particle), intent(in) :: particles(:)
    character(kind=CK, len=3), allocatable :: descr(:)
    type(json_value), pointer :: coordinates_json
    type(json_value), pointer :: particle_json
    integer :: i
    character(kind=CK, len=:), allocatable :: str
    if (size(particles) == 0) return
    call particles(1)%typestr(str)
    call json_add(json_val, "type", str)
    call json_add(json_val, "n", size(particles))
    call particles(1)%description(descr)
    call json_add(json_val, "description", descr)
    call json_create_array(coordinates_json, 'coordinates')
    do i = 1, size(particles)
       call particles(i)%coordinates_to_json(particle_json)
       call json_add(coordinates_json, particle_json)
    end do
    call json_add(json_val, coordinates_json)
  end subroutine particlearray_to_json
  

  !> Performs a trial MC move of @p this particle.
  !!
  !! @param genstates are the random number generator states.
  !! @param simbox is the simulation box in which the @p particle
  !!        resides.
  !! @param temperature is the simulation temperature.
  !! @param nbrs is the neighbours of @p this.
  !! @param pair_ias is the array of interactions corresponding to
  !!        @p nbrs.
  !! @param single_ia is the single-particle interaction acting on
  !!        @p this.
  !! @param dE the change of potential energy for @p this particle due
  !!        to the move.
  !! @param n_trials is the number of trial moves attempted.
  !! @param n_accepted the number of trial moves accepted.
  !! 
  subroutine moveparticle_2(this, genstates, simbox, temperature, nbrs, &
     pair_ias, single_ia, dE, n_trials, n_accepted)
    class(particle), intent(inout) :: this
    class(particlearray_wrapper), intent(in) :: nbrs(:) 
    type(poly_box), intent(in) :: simbox
    real(dp), intent(in) :: temperature
    type(rngstate), intent(inout) :: genstates(0:)
    type(pair_interaction_ptr), intent(in) :: pair_ias(:)
    class(single_interaction), pointer, intent(in) :: single_ia
    real(dp), intent(out) :: dE
    integer, intent(out) :: n_trials, n_accepted
    
    class(particle), allocatable :: newparticle
    class(particle), allocatable :: oldparticle
    integer :: err
    real(dp) :: enew
    real(dp) :: eold
    logical :: isaccepted
    
    enew = 0._dp
    eold = 0._dp
    dE = 0._dp
    isaccepted = .false.
    allocate(newparticle, source=this)
    call newparticle%move(genstates(0))
    call newparticle%set_position(minimage(simbox, newparticle%position()))
    allocate(oldparticle, source=this)
    call this%downcast_assign(newparticle)
    
    call this%energy(nbrs, pair_ias, simbox, single_ia, enew, err)
    
    call this%downcast_assign(oldparticle)
    if(err == 0) then 
       call this%energy(nbrs, pair_ias, simbox, single_ia, &
            enew, err)
       call acceptchange(eold, enew, temperature, genstates(0), isaccepted)
       if(isaccepted) then
          call this%downcast_assign(newparticle)
          dE = enew - eold
       end if
    end if
    
    if (isaccepted) then
       n_accepted = 1
    else
       n_accepted = 0
    end if
    n_trials = 1
  end subroutine moveparticle_2


  !> Returns the potential @p energy of @p this particle due to pairwise
  !! interactions @p pair_ias with @p nbrs and external field
  !! @p single_ia. @p simbox is the simulation box in which @p this
  !! particle and @p nbrs reside. err /= 0 if an error or overlap occurs
  !! in any of the interactions.
  subroutine particle_potential(this, nbrs, pair_ias, simbox, &
       single_ia, energy, err)
    class(particle), intent(in) :: this
    class(particlearray_wrapper), intent(in) :: nbrs(:)
    type(pair_interaction_ptr), intent(in) :: pair_ias(:)
    type(poly_box), intent(in) :: simbox
    class(single_interaction), pointer, intent(in) :: single_ia
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    real(dp) :: e
    integer :: i
    real(dp) :: rcutoff
    real(dp), allocatable :: rijs(:, :)
    logical, allocatable :: cutoff_mask(:)
    integer :: i_nbr
    logical :: use_min_image
    energy = 0
    err = 0

    do i_nbr = 1, size(nbrs)
       use_min_image = .false.
       rcutoff = pair_ias(i_nbr)%ptr%get_cutoff()
       ! Minimum image can only be needed if the cutoff radius goes
       ! outside the box in a periodic dimension.
       if (simbox%xperiodic) then
          if (abs(-0.5 * simbox%lx - this%x) < rcutoff .or. &
               abs(0.5 * simbox%lx - this%x) < rcutoff) use_min_image = .true.
       end if
       if (simbox%yperiodic) then
          if (abs(-0.5 * simbox%ly - this%y) < rcutoff .or. &
               abs(0.5 * simbox%ly - this%y) < rcutoff) use_min_image = .true.
       end if
       if (simbox%zperiodic) then
          if (abs(-0.5 * simbox%lz - this%z) < rcutoff .or. &
               abs(0.5 * simbox%lz - this%z) < rcutoff) use_min_image = .true.
       end if
             
       
       if (err /= 0) exit
       ! Select the calculated interactions by cutoff already here
       if (allocated(rijs)) deallocate(rijs)
       allocate(rijs(3, nbrs(i_nbr)%n))
       if (use_min_image) then
          do i = 1, nbrs(i_nbr)%n
             rijs(:, i) = minimage(simbox, nbrs(i_nbr)%arr(i)%x-this%x,& 
                  nbrs(i_nbr)%arr(i)%y-this%y, nbrs(i_nbr)%arr(i)%z-this%z)     
          end do !! ifort does not vectorize
       else
          do i = 1, nbrs(i_nbr)%n
             rijs(:, i) = [nbrs(i_nbr)%arr(i)%x-this%x,& 
                  nbrs(i_nbr)%arr(i)%y-this%y, nbrs(i_nbr)%arr(i)%z-this%z]
          end do !! ifort does not vectorize
       end if
       if (allocated(cutoff_mask)) deallocate(cutoff_mask)
       allocate(cutoff_mask(nbrs(i_nbr)%n))
       do i = 1, nbrs(i_nbr)%n
          cutoff_mask(i) = dot_product(rijs(:, i), rijs(:, i)) < rcutoff**2
       end do
       cutoff_mask = cutoff_mask .and. nbrs(i_nbr)%mask
       
       do i = 1, nbrs(i_nbr)%n
          if (cutoff_mask(i)) then
             call pair_ias(i_nbr)%ptr%pair_potential(this, nbrs(i_nbr)%arr(i), &
                  rijs(:, i), e, err)
             if (err /= 0) then
                energy = huge(energy)
                return
             else
                energy = energy + e
             end if
          end if
       end do
    end do
    if (err == 0 .and. associated(single_ia)) then
       call single_ia%potential(this, simbox, e, err)
       if (err == 0) energy = energy + e
    end if
  end subroutine particle_potential
  
  
  !> Returns the position of @p aparticle as a vector.
  pure function position(this)
    class(particle), intent(in) :: this
    real(dp), dimension(3) :: position
    position = [this%x, this%y, this%z]
  end function position

  
  !> Assigns the position @p vec to @p aparticle.
  pure subroutine particle_set_position(this, vec)
    class(particle), intent(inout) :: this
    real(dp), dimension(3), intent(in) :: vec
    this%x = vec(1)
    this%y = vec(2)
    this%z = vec(3)
  end subroutine particle_set_position
  

  !> Finalizes @p this.
  subroutine wrapper_delete(this)
    type(particlearray_wrapper), intent(inout) :: this
    if (allocated(this%arr)) deallocate(this%arr)
    if (allocated(this%mask)) deallocate(this%mask)
    this%n = 0
  end subroutine wrapper_delete


  !> Manually finalizes @p this.
  subroutine wrapper_manual_delete(this)
    class(particlearray_wrapper), intent(inout) :: this
    call wrapper_delete(this)
  end subroutine wrapper_manual_delete
  

  !> Assigns @p src to @p this.
  subroutine wrapper_assign(this, src)
    class(particlearray_wrapper), intent(inout) :: this
    type(particlearray_wrapper), intent(in) :: src
    if(.not. allocated(src%arr)) stop &
         'ERROR: wrapper_assign: src%arr not allocated!'
    if (.not. allocated(src%mask)) stop &
         'ERROR: wrapper_assign: src%mask not allocated!'
    ! @todo This could possibly be optimized by not making the (re)allocation
    !       every time, but only when size increases or type of particles
    !       changes.
    if (allocated(this%arr)) deallocate(this%arr)
    allocate(this%arr(size(src%arr)), source=src%arr)
    if (allocated(this%mask)) deallocate(this%mask)
    allocate(this%mask(size(src%mask)), source=src%mask)
    this%n = src%n
  end subroutine wrapper_assign
  
end module m_particle
