!> Implements the domain decomposition algorithm for moving particles.
!! The algorithm relies on the cell list implementation in
!! class_simplelist.
module m_nvt_engine
  use num_kind, only: dp
  use iso_fortran_env, only: output_unit, error_unit
  use class_poly_box, only: poly_box, getx, gety, getz, minimage, &
       isxperiodic, isyperiodic, iszperiodic
  use m_point, only: point, &
       pair_interaction, pair_interaction_ptr, &
       single_interaction, single_interaction_ptr, &
       particlearray_wrapper, particlearray_to_json
  use utils, only: splitstr, join, acceptchange
  !$ use omp_lib
  use class_simplelist, only: simplelist, new_simplelist, simplelist_update, &
       simplelist_nbr_cells, flat_index, simplelist_deallocate, &
       simplelist_nbrmask, simplelist_cell_nbrmask
  use mt_stream, only: rng=>genrand_double1_s, rngstate=>mt_state
  use json_module
  use m_particlejson_parser, only: particlearray_from_json
  use m_json_wrapper, only: get_parameter
  use particle_mover, only: setmaxmoves, getmaxmoves, get_max_translation
  use m_interaction_factory, only: create_pair_interaction, &
       create_single_interaction
  implicit none  

  !> The particlegroup to simplify handling multiple cell lists and and
  !! arrays of particles.
  type particlegroup
     !> The name of the group
     character(len=:), allocatable :: name
     !> The particles in this group.
     class(point), allocatable :: particles(:)
     !> The cell list
     type(simplelist) :: sl
   contains
     procedure :: to_json => particlegroup_to_json
     procedure :: scalepositions
     procedure :: check_particles => particlegroup_check_particles
     final :: particlegroup_finalize
  end type particlegroup

  !> Wrapper to be used for arrays of particlegroups of different types.
  type particlegroup_ptr
     class(particlegroup), pointer :: ptr => null()
  end type particlegroup_ptr

  !> Constructor interface.
  interface particlegroup
     module procedure create_particlegroup
  end interface particlegroup


  !> Stores particles in one cell and its neighbour cells.
  type, extends(particlearray_wrapper) :: domain
     integer :: n_cell = 0
       !! Number of particles in the cell.
  contains
    procedure :: domain_assign
    generic :: assignment(=) => domain_assign
    procedure :: delete => domain_manual_delete
    final :: domain_delete
  end type

  !> Counter for trial particle moves.
  integer, save :: nmovetrials = 0
  
  !> The number of accepted particle moves
  integer, save :: nacceptedmoves = 0
  
  !> The desired acceptance ratio for trial particle moves.
  real(dp), save :: moveratio
  
  !> The simulation temperature.
  real(dp), save :: temperature = -1._dp

  !> Stores the names of the groups that are used in the simulation.
  character(len=20), allocatable, save :: group_names(:)

  !> Stores the particle groups.
  type(particlegroup_ptr), allocatable, save :: groups(:)

  !> The simulation box, which contains the particle groups.
  type(poly_box), save :: simbox
  
  !> Stores a matrix of pair interactions, corresponding to group_names.
  !! if particle groups with names group_names(k) and group_names(l) interact
  !! with each other, then pair_interactions(k, l)%ptr and 
  !! pair_interactions(l, k)%ptr contain that interaction.
  type(pair_interaction_ptr), allocatable, save :: pair_interactions(:, :)

  !> Stores an array of single particle interactions, corresponding to
  !! group_names. 
  type(single_interaction_ptr), allocatable, save :: single_interactions(:)

contains

  !> Initializes the module by reading variables from JSON. Must be
  !! called before any other routine in this module.
  !!
  !! @param json_val contains the JSON.
  !! 
  subroutine nvt_engine_init_json(json_val)
    type(json_value), intent(in), pointer :: json_val
    type(json_value), pointer :: box_json
    type(json_value), pointer :: groups_json => null(), groups_element => null()
    character(len=:, kind=CK), allocatable :: group_name
    integer :: i
    call get_parameter(json_val, 'move_ratio', moveratio, error_lb=0._dp, &
         error_ub=1._dp)
    call get_parameter(json_val, 'temperature', temperature, error_lb=0._dp)
    call get_parameter(json_val, 'nmovetrials', nmovetrials, error_lb=0)
    call get_parameter(json_val, 'nacceptedmoves', nacceptedmoves, &
         error_ub=nmovetrials, error_lb=0)

    ! Getting group names already here, so that we can set up interactions.
    call json_get(json_val, 'particle_groups', groups_json)
    allocate(group_names(json_count(groups_json)))
    group_names = ''
    do i = 1, json_count(groups_json)
       call json_get_child(groups_json, i, groups_element)
       group_name = ''
       call get_parameter(groups_element, 'name', group_name)
       if (group_name == '') then
          write(error_unit, *) 'ERROR: particle_groups(', i, ') has no name!'
          write(error_unit, *) 'Add name to the particle_group!'
          stop
       end if
       if (i > 1) then
          if (any(group_name == group_names(:i-1))) then
             write(error_unit, *) 'Two groups must not have the same name: ' //&
                  group_name
             stop
          end if
       end if
       group_names(i) = group_name
    end do
  
    call create_pair_interactions_json(json_val, group_names, &
         pair_interactions)
    call create_single_interactions_json(json_val, group_names, &
         single_interactions) 

    !! Create box
    call json_get(json_val, 'box', box_json)
    call simbox%from_json(box_json)

    call create_groups(json_val, simbox, pair_interactions, groups)
  end subroutine nvt_engine_init_json


  
  !> Serializes the module state to JSON. The JSON is added to
  !! @p json_val.
  !! 
  subroutine nvt_engine_to_json(json_val)
    type(json_value), pointer, intent(inout) :: json_val
    type(json_value), pointer :: pair_ia_json, pair_ia_element, group_json, &
         group_json_element, box_json, single_ia_json, single_ia_element
    integer :: i, j
    call json_add(json_val, 'temperature', temperature)
    call json_add(json_val, 'move_ratio', moveratio)
    call json_add(json_val, 'nmovetrials', nmovetrials)
    call json_add(json_val, 'nacceptedmoves', nacceptedmoves)
    if (nmovetrials > 0) then
       call json_add(json_val, 'current_move_ratio', &
            real(nacceptedmoves, dp)/real(nmovetrials, dp))
    else
       call json_add(json_val, 'current_move_ratio', 'nan')
    end if
    !! Write pair interactions
    call json_create_array(pair_ia_json, 'pair_interactions')
    do j = 1, size(pair_interactions, 2)
       do i = 1, j
          call json_create_object(pair_ia_element, '')
          call json_add(pair_ia_element, 'participants', &
               [group_names(i), group_names(j)])
          call pair_interactions(i, j)%ptr%to_json(pair_ia_element)
          call json_add(pair_ia_json, pair_ia_element)
       end do
    end do
    call json_add(json_val, pair_ia_json) 
    
    !! Write single interactions
    call json_create_array(single_ia_json, 'single_interactions')
    do j = 1, size(single_interactions)
       if (associated(single_interactions(j)%ptr)) then
          call json_create_object(single_ia_element, '')
          call json_add(single_ia_element, 'participant', group_names(j))
          call single_interactions(j)%ptr%to_json(single_ia_element)
          call json_add(single_ia_json, single_ia_element)
       end if
    end do
    if (json_count(single_ia_json) > 0) then
       call json_add(json_val, single_ia_json) 
    end if
    
    !! Write box
    call json_create_object(box_json, 'box')
    call simbox%to_json(box_json)
    call json_add(json_val, box_json)
    
    !! Write groups
    call json_create_array(group_json, 'particle_groups')
    do i = 1, size(groups)
       call json_create_object(group_json_element, '')
       call groups(i)%ptr%to_json(group_json_element)
       call json_add(group_json, group_json_element)
    end do
    call json_add(json_val, group_json)
    
  end subroutine nvt_engine_to_json

  subroutine nvt_engine_finalize()
    integer :: i, j
    !! Delete interactions
    do i = 1, size(single_interactions)
       if (associated(single_interactions(i)%ptr)) then
          deallocate(single_interactions(i)%ptr)
          single_interactions(i)%ptr => null()
       end if
    end do
    deallocate(single_interactions)
    do j = 1, size(pair_interactions, 2)
       do i = 1, size(pair_interactions, 1)
          if (associated(pair_interactions(i, j)%ptr)) then
             deallocate(pair_interactions(i, j)%ptr)
             pair_interactions(i, j)%ptr => null()
             pair_interactions(j, i)%ptr => null()
          end if
       end do
    end do
    deallocate(pair_interactions)
    !! Delete groups
    do i = 1, size(groups)
       if (associated(groups(i)%ptr)) then
          deallocate(groups(i)%ptr)
          groups(i)%ptr => null()
       end if
    end do
    deallocate(groups)
  end subroutine nvt_engine_finalize

  subroutine nvt_engine_update_decomposition
    integer :: i
    do i = 1, size(groups)
       groups(i)%ptr%sl%min_length = &
            pair_interactions(1, 1)%ptr%get_cutoff() + 2._dp * &
            get_max_translation()
    end do
  end subroutine nvt_engine_update_decomposition
  

  !> Check that @p simbox is large enough if it is periodic.
  subroutine check_simbox(simbox)
    type(poly_box), intent(in) :: simbox
    real(dp) :: largest_cutoff
    integer :: i, j
    largest_cutoff = 0._dp
    do j = 1, size(pair_interactions, 2)
       do i = 1, size(pair_interactions, 1)
          largest_cutoff = max(largest_cutoff, &
               pair_interactions(i, j)%ptr%get_cutoff())
       end do
    end do
    if (simbox%xperiodic .and. getx(simbox) < 2._dp * largest_cutoff) &
         stop 'Simulation box too small!'
    if (simbox%yperiodic .and. gety(simbox) < 2._dp * largest_cutoff) &
         stop 'Simulation box too small!'
    if (simbox%zperiodic .and. getz(simbox) < 2._dp * largest_cutoff) &
         stop 'Simulation box too small!'
  end subroutine check_simbox



  
  !> Reads and creates the particlegroups from the JSON in json_val.
  !!
  !! @param json_val the json_value containing the list of 
  !!        particlegroups.
  !! @param simbox the simulation box.
  !! @param group_names the names of the groups to be read from JSON.
  !! @param pair_interactions the matrix of interactions to be connected
  !!        to the groups.
  !! @param groups the particlegroups read from JSON.
  !!
  subroutine create_groups(json_val, simbox, pair_interactions, &
       groups)
    type(json_value), pointer, intent(in) :: json_val
    type(poly_box), intent(in) :: simbox
    character(kind=CK, len=:), allocatable :: group_name
    type(pair_interaction_ptr), intent(in) :: pair_interactions(:, :)
    type(particlegroup_ptr), allocatable, intent(out) :: groups(:)
    integer :: i, j
    real(dp) :: max_cutoff = 0._dp
    type(json_value), pointer :: groups_json, groups_json_element
    class(point), allocatable :: particles(:)
    !! Find maximum cutoff radius for pair interactions to be used when 
    !! determining the cell sizes of the cell list.
    do i = 1, size(pair_interactions, 2)
       do j = 1, size(pair_interactions, 1)
          max_cutoff = max(max_cutoff, &
          pair_interactions(j, i)%ptr%get_cutoff())
       end do
    end do
    
    call json_get(json_val, 'particle_groups', groups_json)
    allocate(groups(json_count(groups_json)))
    do i = 1, json_count(groups_json)
       call json_get_child(groups_json, i, groups_json_element)
       group_name = ''
       call get_parameter(groups_json_element, 'name', group_name)
       if (group_name == '') then
          write(error_unit, *) 'ERROR: particle_groups(', i, ') has no name!'
          write(error_unit, *) 'Add name to the particle_group!'
          stop
       end if
       call particlearray_from_json(groups_json_element, particles)
       allocate(groups(i)%ptr, source=particlegroup(simbox, particles, &
            min_cell_length=max_cutoff + 2 * get_max_translation(), &
            min_boundary_width=2 * get_max_translation(), name=group_name))
       if (allocated(particles)) deallocate(particles)
    end do
  end subroutine create_groups
  
  
  !> Deserialize pair interactions from JSON. The result is a symmetric
  !! matrix of pair_interaction_ptr:s in which the element pair_ias(i, j)
  !! is the interaction between @p group_names(i) and @p group_names(j).
  !!
  !! @param json_val contains the JSON.
  !! @param group_names the names of the groups used in the simulation.
  !! @param pair_ias the matrix of pair interactions in which the indices
  !!        correspond to group_names. 
  subroutine create_pair_interactions_json(json_val, group_names, pair_ias)
    type(json_value), pointer, intent(in) :: json_val
    character(len=*), intent(in) :: group_names(:)
    type(pair_interaction_ptr), intent(inout), allocatable :: pair_ias(:, :)
    type(json_value), pointer :: pair_ia_json, pair_ia_element
    character(len=len(group_names)), allocatable :: participants(:)
    logical :: found
    integer :: i, j, k, l
    call json_get(json_val, 'pair_interactions', pair_ia_json, found)
    if (found) then
       allocate(pair_ias(size(group_names), size(group_names)))
       do i = 1, json_count(pair_ia_json)
          call json_get_child(pair_ia_json, i, pair_ia_element)
          call get_parameter(pair_ia_element, "participants", participants)
          if (size(participants) == 2) then
             k = 0
             l = 0
             do j = 1, size(group_names)
                if (group_names(j) == participants(1)) k = j
                if (group_names(j) == participants(2)) l = j
             end do
             if (k > 0 .and. l > 0) then
                if (associated(pair_ias(k, l)%ptr)) then
                   write(error_unit, *) &
                        'Warning: only the first interaction with participants ',&
                        participants, ' is used.'
                else
                   pair_ias(k, l)%ptr => create_pair_interaction(pair_ia_element)
                   if (k /= l) then
                      pair_ias(l, k)%ptr => pair_ias(k, l)%ptr
                   end if
                end if
             else
                write(error_unit, *) &
                     'ERROR: Pair interaction has invalid participant.'
                stop 'create_pair_interactions_json unable to continue.'
             end if
          else
             write(error_unit, *) 'ERROR: Pair interaction has ', &
                  size(participants), ' participants.'
             stop 'create_pair_interactions_json unable to continue.'
          end if
       end do
    else
       write(error_unit, *) 'ERROR: pair_interactions not found in json:'
       call json_print(json_val, error_unit)
       stop
    end if
    do j = 1, size(group_names)
       do i = 1, size(group_names)
          if (.not. associated(pair_ias(i, j)%ptr)) then
             write(error_unit, *) 'Interaction between ', group_names(i), &
                  " and ", group_names(j), 'is not defined.'
             stop 'create_pair_interactions_json unable to continue.'
          end if
       end do
    end do
  end subroutine create_pair_interactions_json
  
  
  !> Deserializes single_interactions from JSON. At return, single_ias is
  !! an array, where @p single_ias(i) is the single_interaction
  !! concerning the group named @p group_names(i).
  !!
  !! @param json_val contains the JSON.
  !! @param group_names the names of the groups used in the simulation.
  !! @param single_ias the array of single_interactions.
  !!
  subroutine create_single_interactions_json(json_val, group_names, single_ias)
    type(json_value), pointer, intent(in) :: json_val
    character(len=*), intent(in) :: group_names(:)
    type(single_interaction_ptr), intent(inout), allocatable :: single_ias(:)
    type(json_value), pointer :: single_ia_json, single_ia_element
    character(kind=CK, len=:), allocatable :: participant
    logical :: found
    integer :: i, j, k
    call json_get(json_val, 'single_interactions', single_ia_json, found)
    allocate(single_ias(size(group_names)))
    if (found) then
       do i = 1, json_count(single_ia_json)
          call json_get_child(single_ia_json, i, single_ia_element)
          call get_parameter(single_ia_element, "participant", participant)
          k = 0
          do j = 1, size(group_names)
             if (group_names(j) == participant) k = j
          end do
          if (k > 0) then
             if (associated(single_ias(k)%ptr)) then
                write(error_unit, *) &
                     'Warning: only the first single interaction with ' // &
                     'participant ', participant, ' is used.'
             else
                single_ias(k)%ptr => create_single_interaction(single_ia_element)
             end if
          else
             write(error_unit, *) &
                  'ERROR: Single interaction has invalid participant.'
             stop 'create_single_interactions_json unable to continue.'
          end if
       end do
    end if
  end subroutine create_single_interactions_json
  
  
  !> Assignment operator implementation.
  subroutine domain_assign(this, src)
    class(domain), intent(inout) :: this
    type(domain), intent(in) :: src
    this%particlearray_wrapper = src%particlearray_wrapper
    this%n_cell = src%n_cell
  end subroutine domain_assign

  !> Final procedure.
  subroutine domain_delete(this)
    type(domain), intent(inout) :: this
    this%n_cell = 0
  end subroutine domain_delete

  !> Manual destructor when the final procedure is not applicable.
  subroutine domain_manual_delete(this)
    class(domain), intent(inout) :: this
    call this%particlearray_wrapper%delete()
    this%n_cell = 0
  end subroutine domain_manual_delete

  !> A wrapper around make_particle_moves.
  !!
  !! @param genstates the random number generator states for each
  !!        thread.
  !! @param dE the change in energy after all trial moves have been
  !!        completed.
  !!
  subroutine nvt_update(genstates, dE)
    type(rngstate), intent(inout) :: genstates(0:)
    real(dp), intent(out) :: dE
    integer :: n_trials, n_accepted
    call make_particle_moves(genstates, temperature, dE, n_trials, n_accepted)
    nmovetrials = nmovetrials + n_trials
    nacceptedmoves = nacceptedmoves + n_accepted
  end subroutine nvt_update

  !> Schedules parallel moves of particles using OpenMP with a domain 
  !! decomposition algorithm. 
  !!
  !! @see e.g. G. Heffelfinger and M. Lewitt. J. Comp. Chem., 17(2):250–265,
  !! 1996.
  !!
  !! @param genstates the random number generator states for each
  !!        thread.
  !! @param dE the change in energy after all trial moves have been
  !!        completed.
  !! @param n_trials the number of trial moves attempted.
  !! @param n_accepted the number of accepted moves.
  !! 
  subroutine make_particle_moves(genstates, temperature, dE, n_trials, &
       n_accepted)
    type(rngstate), intent(inout) :: genstates(0:)
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: dE
    integer, intent(out) :: n_accepted, n_trials
    !$ integer :: n_threads
    integer :: thread_id
    !! You may be tempted to make the allocatable arrays automatic, but there's
    !! no performance gained and depending on the system large automatic arrays
    !! may cause a stack overflow.
    real(dp) :: dE_d
    integer :: ix, iy, iz, jx, jy, jz, n_trials_d, n_accepted_d, i_group, i,&
         n_max
    type(domain), allocatable :: ds(:)
    do i_group = 1, size(groups)
#ifdef DEBUG
       if (.not. groups(i_group)%ptr%check_particles(simbox)) then
          write(error_unit, *) 'groups(', i_group, &
               '): some particles are not in the box before moves!'
          stop
       end if
#endif
       call simplelist_update(groups(i_group)%ptr%sl, simbox, &
            groups(i_group)%ptr%particles)
    end do
    thread_id = 0
    dE = 0._dp
    n_accepted = 0
    n_trials = 0
    if (size(groups) == 0) return
    !$ n_threads = 1
    !! Loop over cells. This can be thought of as looping through a 
    !! 2 x 2 x 2 cube of cells.
    !$OMP PARALLEL shared(groups, simbox, genstates, pair_interactions, single_interactions)& 
    !$OMP& private(thread_id, n_threads, n_trials_d, n_accepted_d, ds, dE_d, &
    !$OMP& i_group, n_max)&
    !$OMP& reduction(+:dE, n_accepted, n_trials) 
    !$ thread_id = omp_get_thread_num()
    allocate(ds(size(groups)))
    do i_group = 1, size(groups)
      n_max = min(maxval(groups(i_group)%ptr%sl%counts) * 27, &
        size(groups(i_group)%ptr%particles))
      allocate(ds(i_group)%arr(n_max), &
        source=groups(i_group)%ptr%particles)
      allocate(ds(i_group)%mask(n_max), source = .false.)
      ds(i_group)%n_cell = 0
      ds(i_group)%n = 0
    end do 
    do iz=0, min(1, groups(1)%ptr%sl%nz-1)
       do iy=0, min(1, groups(1)%ptr%sl%ny-1)
          do ix=0, min(1, groups(1)%ptr%sl%nx-1)
             !$OMP DO collapse(3) schedule(dynamic)
             do jz = iz, groups(1)%ptr%sl%nz - 1, 2
                do jy = iy, groups(1)%ptr%sl%ny - 1, 2
                   do jx = ix, groups(1)%ptr%sl%nx - 1, 2
                      
                      !! Collect temp_particles for all groups
                      do i_group = 1, size(groups)
                         call set_domain(groups(i_group)%ptr, simbox, &
                              jx, jy, jz, ds(i_group))
                      end do
                      !! Move particles
                      do i_group = 1, size(groups)
                         call domain_move(ds, i_group, &
                              genstates(thread_id:thread_id), simbox, &
                              temperature, pair_interactions, &
                              single_interactions, dE_d, n_trials=n_trials_d, &
                              n_accepted=n_accepted_d)
                         dE = dE + dE_d
                         n_trials = n_trials + n_trials_d
                         n_accepted = n_accepted + n_accepted_d
                      end do
                      !! Synchronize
                      !! :TODO: Generalize to a sync domain operation?
                      do i_group = 1, size(groups)
                         do i = 1, ds(i_group)%n_cell
                            call groups(i_group)%ptr%particles(&
                                 groups(i_group)%ptr%sl%indices(i, jx, jy, jz)&
                                 )%downcast_assign(ds(i_group)%arr(i))
                         end do
                      end do
                   end do
                end do
             end do
             !$OMP END DO 
             !! The end of parallelized loop forces an implicit barrier.
             !! Memory view of the threads is also synchronized here.
          end do
          !$OMP BARRIER
       end do
       !$OMP BARRIER
    end do
    do i_group = 1, size(groups)
       call ds(i_group)%delete()
    end do   
    !$OMP END PARALLEL

    do i_group = 1, size(groups)
#ifdef DEBUG
       if (.not. groups(i_group)%ptr%check_particles(simbox)) then
          write(error_unit, *) 'groups(', i_group, &
               '): some particles are not in the box!'
          stop
       end if
#endif
       call simplelist_update(groups(i_group)%ptr%sl, simbox, &
            groups(i_group)%ptr%particles)
    end do
  end subroutine make_particle_moves
  
  !> Resets the counters used when adjusting moves.
  subroutine nvt_engine_reset_counters()
    nacceptedmoves = 0
    nmovetrials = 0
  end subroutine nvt_engine_reset_counters

  !> Assigns the cell @p jx @p jy @p jz to the domain @p d.
  !!
  !! @param this the particlegroup of the domain.
  !! @param simbox the simulation box.
  !! @param jx, jy, jz the cell index.
  !! @param d the domain with particles from cell jx, jy, jz at return.
  !!
  subroutine set_domain(this, simbox, jx, jy, jz, d)
    type(particlegroup), intent(in) :: this
    type(poly_box), intent(in) :: simbox
    integer, intent(in) :: jx, jy, jz
    type(domain), intent(inout) :: d
    logical, allocatable :: nbr_mask(:)
    integer :: i
    integer, allocatable :: helper(:)
    d%mask = .false.
    allocate(nbr_mask(size(this%particles)), source=.false.)
    d%n_cell = this%sl%counts(jx, jy, jz)
    call simplelist_cell_nbrmask(this%sl, simbox, jx, jy, jz, nbr_mask)
    d%n = count(nbr_mask)
    !! Remove particles in jx, jy, jz from mask:
    nbr_mask(this%sl%indices(1:d%n_cell, jx, jy, jz)) = .false.
    !! GFortran 6.0.0 is not fine with pack from a polymorphic
    !! array class(point) so we'll use a workaround.
    allocate(helper, source=[this%sl%indices(1:d%n_cell, jx, jy, jz), &
         pack([(i, i = 1, size(this%particles))], nbr_mask)])
    if (allocated(nbr_mask)) deallocate(nbr_mask)
    do i = 1, d%n
       call d%arr(i)%downcast_assign(this%particles(helper(i)))
    end do
    d%mask(1:d%n) = .true.
  end subroutine set_domain  

  
  !> Moves particle in the cell of the domain @p domains(i_d). Used in
  !! the domain decomposition algorithm.
  !!
  !! @param domains are the subsystems of particles corresponding to
  !! 	    groups given to the make_particle_moves algorithm.
  !! @param i_d is the index of the domain to be moved in @p domains.
  !! @param genstates the random number generator states.
  !! @param simbox the simulation box.
  !! @param temperature the simulation temperature.
  !! @param pair_interactions the pair interactions concerning domains(i_d)
  !! @param single_interactions the single-particle interaction concerning 
  !! 	    domains(i_d). 
  !! @param dE the total change in energy due to the moves.
  !! @param n_trials the number of trial moves that were performed.
  !! @param n_accepted the number of accepted moves of n_trials.
  !!
  subroutine domain_move(domains, i_d, genstates, simbox, temperature, &
       pair_interactions, single_interactions, dE, n_trials, n_accepted)
    type(domain), intent(inout) :: domains(:)
    integer, intent(in) :: i_d
    type(rngstate), intent(inout) :: genstates(:)
    type(poly_box), intent(in) :: simbox
    real(dp), intent(in) :: temperature
    type(pair_interaction_ptr), intent(in) :: pair_interactions(:, :)
    type(single_interaction_ptr) :: single_interactions(:)
    real(dp), intent(out) :: dE
    integer, intent(out), optional :: n_trials, n_accepted
    integer :: j
    class(point), allocatable :: newparticle
    integer :: err
    real(dp) :: enew
    real(dp) :: eold
    logical :: isaccepted
    n_accepted = 0
    n_trials = 0
    dE = 0.
    if (domains(i_d)%n_cell > 0) allocate(newparticle, &
         source=domains(i_d)%arr(1))
    do j = 1, domains(i_d)%n_cell
       domains(i_d)%mask(j) = .false.
#ifdef DEBUG
       call newparticle%downcast_assign(domains(i_d)%arr(j), err)
       if (err /= 0) then
          write(error_unit, *) 'newparticle%downcast_assign: err=', err
          stop 'domain_move unable to continue'
       end if
#else
       call newparticle%downcast_assign(domains(i_d)%arr(j))
#endif
       call newparticle%move(genstates(1))
       call newparticle%set_position(simbox%minimage(newparticle%x, &
            newparticle%y, newparticle%z))
       call newparticle%energy(domains, pair_interactions(:, i_d), simbox, &
            single_interactions(i_d)%ptr, enew, err)
       if(err == 0) then 
          call domains(i_d)%arr(j)%energy(domains, pair_interactions(:, i_d),&
               simbox, single_interactions(i_d)%ptr, eold, err)
          if (err /= 0) then
             write(error_unit, *) 'ERROR: err=', err, ' domain=', i_d, &
                  ' old particle ', j
             if (err == 1) write(error_unit, *) &
                  'Particle energy calculation resulted in overlap before move.'
             stop 'Stopped by domain_move.'
          end if
          call acceptchange(eold, enew, temperature, genstates(1), isaccepted)
          if(isaccepted) then
#ifdef DEBUG
             call domains(i_d)%arr(j)%downcast_assign(newparticle, err)
             if (err /= 0) then
                write(error_unit, *) 'domains(', i_d, &
                     '%arr(', j, ')%downcast_assign: err = ', err
                stop 'domain_move unable to continue'
             end if
#else
             call domains(i_d)%arr(j)%downcast_assign(newparticle)
#endif
             dE = dE + enew - eold
             n_accepted = n_accepted + 1
          end if
       end if
       domains(i_d)%mask(j) = .true.
    end do
    n_trials = n_trials + domains(i_d)%n_cell
  end subroutine domain_move

  !> Adjusts the maximum moves of the particles.
  subroutine nvt_engine_update_max_moves()
    real(dp) :: newdthetamax
    real(dp) :: newdximax
    if (nmovetrials > 0) then
       call getmaxmoves(newdximax, newdthetamax)
       !! Adjust translation
       newdximax = newmaxvalue(nmovetrials, nacceptedmoves, moveratio, &
            newdximax)
       !! Update the minimum cell side length of the cell list because the
       !! maximum translation has changed: 
     
       !! This should adjust rotations < pi/2 to move the particle end as
       !! much as a random translation. 4.4 is the assumed molecule 
       !! length-to-breadth ratio.
       newdthetamax = 2._dp * asin(newdximax/4.4_dp) 
       call setmaxmoves(newdximax, newdthetamax)
    end if
  end subroutine nvt_engine_update_max_moves
  
  !> Returns a new trial move parameter value calculated from the desired
  !! acceptance ratio.
  !! 
  !! @param ntrials the total number of trials of the kind of move in
  !!        question.
  !! @param naccepted number of accepted trials.
  !! @param desiredratio the desired acceptance ratio for trial moves.
  !! @param oldvalue the old value of the parameter setting the maximum
  !!        size for the trial move in question.
  !!
  function newmaxvalue(ntrials, naccepted, desiredratio, oldvalue) &
       result(newvalue)
    integer, intent(in) :: ntrials
    integer, intent(in) :: naccepted
    real(dp), intent(in) :: desiredratio
    real(dp), intent(in) :: oldvalue
    real(dp) :: newvalue
    real(dp), parameter :: multiplier = 1.05_dp
    if (real(naccepted, dp) / real(ntrials, dp) > desiredratio) then
       newvalue = oldvalue * multiplier
    else 
       newvalue = oldvalue / multiplier
    end if
  end function newmaxvalue


  ! particle_group routines:

  !> Creates a particlegroup.
  !!
  !! @param simbox the box in which the particles reside.
  !! @param particles the particles in the group.
  !! @param min_cell_length the minimum cell side length in the cell
  !!        list.
  !! @param min_boundary_width the width of the "boundary" area around
  !!        the cell. When any particle crosses the boundary, the cell
  !!        list should be updated.
  !! 
  !! @return the particlegroup which has a cell list with the given
  !!         specs.
  !!
  function create_particlegroup(simbox, particles, min_cell_length, &
       min_boundary_width, name) result(group)
    type(poly_box), intent(in) :: simbox
    class(point), intent(in) :: particles(:)
    real(dp), intent(in) :: min_cell_length, min_boundary_width
    character(kind=CK, len=*), intent(in) :: name
    type(particlegroup) :: group
    group%name = name
    allocate(group%particles(size(particles)), &
         source=particles(1:size(particles)))
    !$ if (.true.) then
    !$ call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
    !$& min_boundary_width, is_nx_even = isxperiodic(simbox), &
    !$& is_ny_even = isyperiodic(simbox), is_nz_even = iszperiodic(simbox))
    !$ else 
    call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
         min_boundary_width)
    !$ end if
  end function create_particlegroup


  !> Serializes the particlegroup @p this to JSON @p json_val.
  subroutine particlegroup_to_json(this, json_val)
    class(particlegroup), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    if (allocated(this%name)) call json_add(json_val, 'name', this%name)
    !! particlearray
    call particlearray_to_json(json_val, this%particles)
  end subroutine particlegroup_to_json

  
  !> Final routine.
  subroutine particlegroup_finalize(group)
    type(particlegroup), intent(inout) :: group
    call simplelist_deallocate(group%sl)
    if(allocated(group%particles)) deallocate(group%particles)
  end subroutine particlegroup_finalize
  

  !> Computes the total @p energy of the system. Interactions are
  !! computed cell-by-cell and group-by-group.
  !!
  !! @param energy the total energy at return.
  !! @param err is non-zero if an error occurs, such as two particles
  !!        being too close to each other.
  !!
  subroutine total_energy(energy, err)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    integer :: i, ix, iy, iz, i_group
    real(dp) :: energy_j
    integer :: nbr_cells(3, 27), n_nbr_cells
    integer :: j_group
    energy = 0._dp
    err = 0
    if (size(groups) == 0) return
    !$OMP PARALLEL default(shared) reduction(+:energy, err)& 
    !$OMP& private(energy_j, i, nbr_cells, n_nbr_cells)
    !$OMP DO collapse(3) schedule(dynamic)
    do ix = 0, groups(1)%ptr%sl%nx - 1 
       do iy = 0, groups(1)%ptr%sl%ny - 1
          do iz = 0, groups(1)%ptr%sl%nz - 1
             if (err == 0) then
                do i_group = 1, size(groups)
                   ! 1. compute inside ix, iy, iz in i_group
                   call cell_energy(groups(i_group)%ptr, ix, iy, iz, simbox, &
                        pair_interactions(i_group, i_group)%ptr, &
                        single_interactions(i_group)%ptr, energy_j, err)
                   if (err /= 0) exit
                   energy = energy + energy_j
                end do
             end if
             if (err == 0) then
                do i_group = 1, size(groups) - 1
                   !! 2. compute with ix, iy, iz in j_group > i_group
                   do j_group = i_group + 1, size(groups)
                      !write(*, *) "compute with ", ix, iy, iz, " in ",
                      ! j_group, " > ", i_group
                      call cell_pair_energy(groups(i_group)%ptr, ix, iy, iz, &
                           groups(j_group)%ptr, ix, iy, iz, simbox, &
                           pair_interactions(i_group, j_group)%ptr, energy_j, err)
                      if (err /= 0) exit
                      energy = energy + energy_j
                   end do
                   if (err /= 0) exit
                end do
             end if
             if (err == 0) then
                ! 3. for all j_group (including i_group) compute where
                ! ix + nx * iy + nx * ny * iz < jx + jy * nx + jz * nx * ny
                ! and jx, jy, jz is a neighbour of ix, iy, iz.
                call simplelist_nbr_cells(groups(i_group)%ptr%sl, ix, iy, iz, &
                     nbr_cells, n_nbr_cells)
                do i_group = 1, size(groups)
                   do i = 1, n_nbr_cells
                      if (flat_index(groups(i_group)%ptr%sl, nbr_cells(1, i), &
                           nbr_cells(2, i), nbr_cells(3, i)) > &
                           flat_index(groups(i_group)%ptr%sl, ix, iy, iz)) then
                         do j_group = 1, size(groups)
                            call cell_pair_energy(groups(i_group)%ptr, &
                                 ix, iy, iz, &
                                 groups(j_group)%ptr, nbr_cells(1, i), &
                                 nbr_cells(2, i), nbr_cells(3, i), simbox, &
                                 pair_interactions(i_group, j_group)%ptr, &
                                 energy_j, err)
                            if (err /= 0) exit
                            energy = energy + energy_j
                         end do
                      end if
                      if (err /= 0) exit
                   end do
                   if (err /= 0) exit
                end do
             end if
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL  
  end subroutine total_energy

  !> The "internal" energy of cell @p ix, @p iy, @p iz in @p this group.
  !!
  !! @param this the particlegroup
  !! @param ix, iy, iz the cell indices.
  !! @param simbox the simulation box.
  !! @param pair_ia the pair interaction between particles in this
  !!        group.
  !! @param single_ia the single interaction concerning this group.
  !! @param energy the energy at return.
  !! @param err is non-zero if any error occurs.
  !!
  subroutine cell_energy(this, ix, iy, iz, simbox, pair_ia, single_ia, &
       energy, err)
    type(particlegroup), intent(in) :: this
    integer, intent(in) :: ix, iy, iz
    type(poly_box), intent(in) :: simbox
    class(pair_interaction), intent(in) :: pair_ia
    class(single_interaction), pointer, intent(in) :: single_ia
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    integer :: i, j
    real(dp) :: energy_ij, rij(3)
    energy = 0
    err = 0
    do i = 1, this%sl%counts(ix, iy, iz) - 1
       associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
         do j = i + 1, this%sl%counts(ix, iy, iz)
            associate(particlej => &
                 this%particles(this%sl%indices(j, ix, iy, iz)))
              rij = minimage(simbox,&
                   particlej%x - particlei%x, particlej%y - particlei%y,&
                   particlej%z - particlei%z)
              if (dot_product(rij, rij) < pair_ia%get_cutoff()**2) then
                 call pair_ia%pair_potential(particlei, particlej, rij, &
                      energy_ij, err)
                 if (err /= 0 ) return
                 energy = energy + energy_ij
              end if
            end associate
         end do
       end associate
    end do
    if (associated(single_ia)) then
       do i = 1, this%sl%counts(ix, iy, iz)
          associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
            call single_ia%potential(particlei, simbox, energy_ij, err)
            if (err /= 0) exit
            energy = energy + energy_ij
          end associate
       end do
    end if
  end subroutine cell_energy

  
  !> The total pair interaction between particles in groups @p this and
  !! @p another in cells [ix, iy, iz] and [jx, jy, jz], respectively.
  !!
  !! @param this the particlegroup
  !! @param ix, iy, iz the cell index for @p this.
  !! @param another the other particlegroup.
  !! @param jx, jy, jz the cell index for @p another.
  !! @param simbox the simulation box in which the particles reside.
  !! @param pair_ia the pair_interaction between @p this and @p another.
  !! @param energy the energy at return.
  !! @param err is non-zero if any pairwise-computed interaction
  !!        results in non-zero err.
  !!
  subroutine cell_pair_energy(this, ix, iy, iz, another, jx, jy, jz, simbox, &
     pair_ia, energy, err)
  type(particlegroup), intent(in) :: this, another
  integer, intent(in) :: ix, iy, iz, jx, jy, jz
  type(poly_box), intent(in) :: simbox
  class(pair_interaction), intent(in) :: pair_ia
  real(dp), intent(out) :: energy
  integer, intent(out) :: err
  integer :: i, j
  real(dp) :: energy_ij, rij(3)
  !write(*, *) "compute ", this%name, ix, iy, iz, " with ", &
  !     another%name, jx, jy, jz
  energy = 0
  if ((ix-jx)**2 + (iy-jy)**2 + (iz-jz)**2 <= 3) then
  ! No need to compute minimum image.
  do i = 1, this%sl%counts(ix, iy, iz)
     associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
       do j = 1, another%sl%counts(jx, jy, jz)
          associate(particlej => another%particles(&
               another%sl%indices(j, jx, jy, jz)))
            rij = [particlej%x - particlei%x, particlej%y - particlei%y, &
                 particlej%z - particlei%z]
            if (dot_product(rij, rij) < pair_ia%get_cutoff()**2) then
               call pair_ia%pair_potential(particlei, particlej, rij, &
                    energy_ij, err)
               if (err /= 0) return
               energy = energy + energy_ij
            end if
          end associate
       end do
     end associate
  end do
  else
  do i = 1, this%sl%counts(ix, iy, iz)
     associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
       do j = 1, another%sl%counts(jx, jy, jz)
          associate(particlej => another%particles(&
               another%sl%indices(j, jx, jy, jz)))
            rij = minimage(simbox, particlej%x - particlei%x, &
                 particlej%y - particlei%y, particlej%z - particlei%z)
            if (dot_product(rij, rij) < pair_ia%get_cutoff()**2) then
               call pair_ia%pair_potential(particlei, particlej, rij, &
                    energy_ij, err)
               if (err /= 0) return
               energy = energy + energy_ij
            end if
          end associate
       end do
     end associate
  end do
  end if
end subroutine cell_pair_energy


!> Scales the positions of @p particles with the same factors that were
!! used to scale the simulation box dimensions from @p oldbox to
!! @p newbox.
!!
!! @param this the particlegroup for which the particles are scaled.
!! @param oldbox, newbox the new and old simulation boxes.
!! 
subroutine scalepositions(this, oldbox, newbox)
  class(particlegroup), intent(inout) :: this
  type(poly_box), intent(in) :: oldbox
  type(poly_box), intent(in) :: newbox
  integer :: i
  !! Don't turn the loop below to array syntax. At least GCC6.0.0 (BETA) can
  !! not handle these arrays correctly.
  do i = 1, size(this%particles)
     this%particles(i)%x = this%particles(i)%x * getx(newbox) / getx(oldbox)
     this%particles(i)%y = this%particles(i)%y * gety(newbox) / gety(oldbox)
     this%particles(i)%z = this%particles(i)%z * getz(newbox) / getz(oldbox)
  end do
  call simplelist_update(this%sl, newbox, this%particles)
end subroutine scalepositions

!> Checks that all particles in @p this group are inside the @p simbox.
!!
!! @param the particlegroup to check.
!! @param simbox the simulation box in which the particles should
!!        reside.
function particlegroup_check_particles(this, simbox) result(res)
  class(particlegroup), intent(in) :: this
  class(poly_box), intent(in) :: simbox
  logical :: res
  integer :: i
  res = .true.
  do i = 1, size(this%particles)
     res = res .and. simbox%check_position(this%particles(i)%position())
  end do
end function particlegroup_check_particles


  
end module m_nvt_engine
