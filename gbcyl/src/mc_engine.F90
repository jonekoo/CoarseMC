!> Implements the steering of the simulation and controlling it's
!! input and output.
module mc_engine
  !$ use omp_lib
  use m_particlegroup, only: total_energy,  particlegroup, particlegroup_ptr
  use m_npt_engine, only: npt_engine_reset_counters, &
       npt_engine_update_max_scaling, npt_engine_finalize, &
       npt_engine_to_json, npt_engine_init_json, pressure, npt_update
  use m_nvt_engine, only: temperature, nvt_update, &
       nvt_engine_update_max_moves, nvt_engine_reset_counters, &
       nvt_engine_to_json, nvt_engine_init_json
  use mt_stream
  use m_fileunit
  use class_poly_box, only: poly_box, getx, gety, getz
  use class_parameterizer
  use class_parameter_writer
  use m_particle, only: particle, pair_interaction_ptr, &
       particlearray_wrapper, single_interaction_ptr
  use particle_mover, only: particlemover_init, particlemover_writeparameters,&
       get_max_translation, getmaxmoves, setmaxmoves, particlemover_to_json
  use beta_exchange, only: write_stats, be_reset_counters, &
       beta_exchange_init => init, be_finalize => finalize, try_beta_exchanges,&
       be_get_acceptance_ratios
  use num_kind, only: dp
  use iso_fortran_env, only: error_unit, output_unit
  use json_module
  use m_json_wrapper, only: get_parameter
  use m_interaction_factory, only: create_pair_interaction, &
       create_single_interaction
  use m_particlejson_parser, only: particlearray_from_json
  implicit none

  private :: mce_init_common
  
  !> The number of equilibration MC sweeps in the simulation.  Equilibration
  !! sweeps are used to do all kinds of adjusting and should be discarded
  !! from the analysis of results.
  integer, save :: nequilibrationsweeps = 0
  
  !> The number of production MC sweeps in the simulation.
  integer, save :: nproductionsweeps = 0
  
  !> The periodicity of writing the configuration of molecules and
  !! simulation parameters on the disk. If productionperiod=400, every
  !! 400th MC sweep the configuration and parameters are written to the
  !! disk. 
  integer, save :: productionperiod = 1
  
  !> The period of adjusting maximum trial move sizes during the
  !! equilibration MC sweeps.
  integer, save :: moveadjustperiod = 100
  
  !> The number of sweeps between parallel tempering updates.
  integer, save :: pt_period = 1
  
  !> If adjusting of temperatures in the temperature series for parallel 
  !! tempering is used defines the period for these actions. (DEPRECATED)
  integer, save :: ptadjustperiod = 100
  
  !> The period of writing a restart point on the disk. This includes the
  !! configuration of molecules, the simulation parameters and the random
  !! number generator state. 
  integer, save :: restartperiod = 10000
  
  !> The sweep counter.
  integer, save :: isweep = 0
  
  !> The input and output unit used for reading and writing the geometry of 
  !! molecules and the simulation box.
  integer, save :: coordinateunit
  
  !> The (MPI) id of this process formatted to a character.
  character(len = 9), save :: idchar

  !> Output file name.
  character(len=:), allocatable, save :: fn_parameters_out

  !> Restart file name.
  character(len=:), allocatable, save :: fn_parameters_restart

  !> The random number generator states.
  type(mt_state), allocatable, save :: mts(:)
  
  !> The random number generator seed.
  integer, save :: seed = -1

  !> Stores the particle groups.
  type(particlegroup_ptr), allocatable, save :: groups(:)

  !> The simulation box, which contains the particle groups.
  type(poly_box), save :: simbox

  !> Stores the names of the groups that are used in the simulation.
  character(len=20), allocatable, save :: group_names(:)

  !> Stores a matrix of pair interactions, corresponding to group_names.
  !! if particle groups with names group_names(k) and group_names(l) interact
  !! with each other, then pair_interactions(k, l)%ptr and 
  !! pair_interactions(l, k)%ptr contain that interaction.
  type(pair_interaction_ptr), allocatable, save :: pair_interactions(:, :)

  !> Stores an array of single particle interactions, corresponding to
  !! group_names. 
  type(single_interaction_ptr), allocatable, save :: single_interactions(:)

  !> The total energy.
  real(dp), save :: etotal = 0._dp

  !> The string defining the ensemble (npt or nvt)
  character(len=:), allocatable, save :: ensemble

  !> The id of this replica.
  integer, save :: replica_id =  0
  integer, save :: n_replicas = 1
  
contains

  subroutine mce_init_json(id, n_tasks, parameter_infile, parameter_outfile, &
       parameter_restartfile)
    !! Initializes the program. Must be called before run or finalize.

    integer, intent(in) :: id
      !! The replica index of the process for parallel tempering.
    integer, intent(in) :: n_tasks
      !! The number of replicas.
    character(len=*), intent(in) :: parameter_infile
      !! The input file containing parameters and coordinates.
    character(len=*), intent(in) :: parameter_outfile
      !! The name of the file to which parameters and coordinates are
      !! logged during the simulation.
    character(len=*), intent(in) :: parameter_restartfile
      !! The name of the file to which the parameters and coordinates
      !! are written for checkpointing/restarting.
    type(json_value), pointer :: json_val, box_json
    logical :: status_ok
    character(len=:), allocatable :: error_msg
    fn_parameters_out = parameter_outfile
    fn_parameters_restart = parameter_restartfile
    write(idchar, '(I9)') id 
    call json_parse(parameter_infile, json_val)
    if (json_failed()) then
       call json_print_error_message(error_unit)
       stop
    end if
    call mce_from_json(json_val)
    if (json_failed()) then
        call json_check_for_errors(status_ok, error_msg)
        write(error_unit, *) 'Error: '//error_msg
        call json_clear_exceptions()
        stop 'json_failed'
    end if
 
    call create_pair_interactions_json(json_val, group_names, &
         pair_interactions)
    call create_single_interactions_json(json_val, group_names, &
         single_interactions)
    call json_get(json_val, 'box', box_json)
    call simbox%from_json(box_json)

    call create_groups(json_val, simbox, pair_interactions, groups)
    
    if (json_failed()) call json_print_error_message(error_unit)
    call mce_init_common(id, n_tasks)
    
  end subroutine mce_init_json


  !> Initializes this module and its dependencies.
  !!
  !! @param id the process id for an MPI process.
  !! @param n_tasks total number of MPI processes.
  !!
  subroutine mce_init_common(id, n_tasks)
    integer, intent(in) :: id
    integer, intent(in) :: n_tasks 
    integer :: err

    replica_id = id
    n_replicas = n_tasks
    call mce_init_rng(id, n_tasks)    
    call beta_exchange_init(1._dp / temperature)
    call total_energy(groups, simbox, pair_interactions, &
         single_interactions, etotal, err)
    if (err /= 0) then
       write(error_unit, *) &
            'ERROR: mce_init: total_energy returned err = ', err
       stop 'Program stopped by mce_init_common.'
    end if
    call makerestartpoint
  end subroutine mce_init_common
  

subroutine mce_from_json(json_val)
  type(json_value), pointer, intent(in) :: json_val
  type(json_value), pointer :: groups_json => null(), groups_element => null()
  character(len=:, kind=CK), allocatable :: group_name
  integer :: i
  call get_parameter(json_val, 'replica_id', replica_id)
  call get_parameter(json_val, 'n_equilibration_sweeps', &
       nequilibrationsweeps, error_lb=0)
  call get_parameter(json_val, 'n_production_sweeps', &
       nproductionsweeps, error_lb=0)
  call get_parameter(json_val, 'production_period', &
       productionperiod, error_lb=1)
  call get_parameter(json_val, 'i_sweep', isweep, error_lb=0)  
  call get_parameter(json_val, 'move_adjusting_period', &
       moveadjustperiod, error_lb=1)
  call get_parameter(json_val, 'pt_period', pt_period, error_lb=1)
  call get_parameter(json_val, 'restartperiod', restartperiod, error_lb=1)
  call get_parameter(json_val, 'seed', seed, warn_lb=1000)
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
  ensemble = 'npt'
  call get_parameter(json_val, 'ensemble', ensemble)
  !! Initialize modules.
  call particlemover_init(json_val)
  call nvt_engine_init_json(json_val)
  if (ensemble == 'npt') then
     call npt_engine_init_json(json_val)
  end if
end subroutine mce_from_json

!> Initializes independen random number generator states for each MPI
!! process and thread.
!!
!! @param id The id of the MPI process.
!! @param n_tasks The total number of MPI processes.
!!
subroutine mce_init_rng(id, n_tasks)
  integer, intent(in) :: id, n_tasks
  integer :: thread_id = 0, n_threads = 1
  type(mt_state) :: temp_mts
  !$ n_threads = omp_get_max_threads()
  
  allocate(mts(0:n_threads - 1))
  call set_mt19937
  !call new(mts(thread_id))
  call new(temp_mts)
  if (seed < 0) then 
     call system_clock(seed)
     write(error_unit, *) "Warning: Seeding RNG with system clock."
  end if
  if (seed < 0) then
     write(error_unit, *) &
          "Warning: system_clock query failed, using default seed 1234567."
     seed = 1234567 
  end if
  call init(temp_mts, seed)
  mts(0) = temp_mts
  
  !$ write(output_unit, *) 'Running with ', n_threads, ' threads.'
  
  !! Give different random number streams to each OpenMP thread inside a
  !! MPI task.
  !$OMP PARALLEL DO
  !$ do thread_id = 1, n_threads-1
  if (id + n_tasks * thread_id > 0) then 
     call create_stream(temp_mts, mts(thread_id), id + n_tasks * thread_id)
  end if
  !$ end do
  !$OMP END PARALLEL DO
end subroutine mce_init_rng


!> Reads and creates the particlegroups from the JSON in json_val.
!!
!! @param json_val the json_value containing the list of particlegroups.
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
  class(particle), allocatable :: particles(:)
  !! Find maximum cutoff radius for pair interactions to be used when 
  !! determining the cell sizes of the cell list.
  do i = 1, size(pair_interactions, 2)
     do j = 1, size(pair_interactions, 1)
        max_cutoff = max(max_cutoff, pair_interactions(j, i)%ptr%get_cutoff())
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

!> Finalizes the simulation.
!!
!! @param id is the MPI process id.
!!
subroutine finalize(id)
  integer, intent(in) :: id
  integer :: i, j
  close(coordinateunit)
  if (mod(isweep, restartperiod) /= 0) then
     ! Make restartpoint at the end only if it was not already made for
     ! this sweep.
     call makerestartpoint
  end if
  if (ensemble == 'npt') then
     call npt_engine_finalize
  end if
  call be_finalize
  do i = 0, size(mts) - 1
     call delete(mts(i))
  end do
  !! Delete interactions
  if (allocated(mts)) deallocate(mts)
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
  if (id == 0) write (output_unit, *) &
       'Program ptgbcyl was finalized succesfully.'
end subroutine 

  
!> Runs one sweep of Metropolis Monte Carlo updates to the system. A
!! full Parallel tempering NPT-ensemble sweep consists of trial moves
!! of particles, trial scaling(s) of the simulation box (barostat) and 
!! trial swaps of temperature with other replicas.
!! 
!! @param simbox the simulation box.
!! @param groups the particlegroups.
!! @param genstates random number generator states for all threads.
!! @param isweep the sweep counter.
!!  
subroutine sweep(simbox, groups, genstates, isweep)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup_ptr), intent(inout) :: groups(:)
  type(mt_state), intent(inout) :: genstates(0:)
  integer, intent(in) :: isweep
  real(dp) :: beta, dE
#ifdef DEBUG
  integer :: err
  write(output_unit, *) 'Move particles'
#endif
  dE = 0.
  call nvt_update(groups, genstates, simbox, pair_interactions, &
       single_interactions, dE)
#ifdef DEBUG
  write(output_unit, *) "etotal + dE = ", etotal, " + ", dE, " = ", etotal + dE
#endif
  etotal = etotal + dE
#ifdef DEBUG
  call total_energy(groups, simbox, pair_interactions, &
       single_interactions, etotal, err)
  write(output_unit, *) "etotal updated = ", etotal
#endif
  if (ensemble == 'npt') then
#ifdef DEBUG
     write(output_unit, *) 'Make volume move(s)'
#endif
     call npt_update(simbox, groups, genstates(0), pair_interactions, &
          single_interactions, etotal)
  end if
  call check_simbox(simbox)
  if (mod(isweep, pt_period) == 0) then
     beta = 1._dp / temperature
     select case (ensemble)
     case ('npt')
        call try_beta_exchanges(beta, etotal + pressure * simbox%volume(), 3,&
             genstates(0)) 
     case ('nvt')
        call try_beta_exchanges(beta, etotal, 3, genstates(0))
     case default
        write(error_unit, *) 'Warning: Unknown ensemble. ' // &
             'Parallel tempering disabled.'
     end select
     temperature = 1._dp / beta
  end if
end subroutine sweep


!> Serializes the module and it's dependencies to JSON value @p
!! json_val.
subroutine mce_to_json(json_val)
  type(json_value), pointer, intent(out) :: json_val
  integer :: i, j
  type(json_value), pointer :: pair_ia_json, pair_ia_element, group_json, &
       group_json_element, box_json, single_ia_json, single_ia_element
  call json_create_object(json_val, 'mc_engine')
  call json_add(json_val, 'replica_id', replica_id)
  call json_add(json_val, 'n_equilibration_sweeps', &
  nequilibrationsweeps)
  call json_add(json_val, 'n_production_sweeps', nproductionsweeps)
  call json_add(json_val, 'i_sweep', isweep)
  call json_add(json_val, 'production_period', productionperiod)
  call json_add(json_val, 'move_adjusting_period', moveadjustperiod)
  call json_add(json_val, 'pt_period', pt_period)
  call json_add(json_val, 'restartperiod', restartperiod)
  call json_add(json_val, 'seed', seed)
  call json_add(json_val, 'volume', simbox%volume())
  call json_add(json_val, 'total_energy', etotal)
  call json_add(json_val, 'enthalpy', etotal + &
         simbox%volume() * pressure)
  call json_add(json_val, 'ensemble', ensemble)

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

  call particlemover_to_json(json_val)
  call nvt_engine_to_json(json_val)
  if (ensemble == 'npt') then
     call npt_engine_to_json(json_val)
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
end subroutine mce_to_json


!> Runs the simulation. 
subroutine run
  integer :: i
  do while (isweep < nequilibrationsweeps + nproductionsweeps)
    isweep = isweep + 1
    call sweep(simbox, groups, mts, isweep)
    if (isweep <= nequilibrationsweeps) then
      if (moveadjustperiod /= 0) then
         if (mod(isweep, moveadjustperiod) == 0) then
            call updatemaxvalues
            do i = 1, size(groups)
               groups(i)%ptr%sl%min_length = &
                    pair_interactions(1, 1)%ptr%get_cutoff() + 2._dp * &
                    get_max_translation()
            end do
         end if
      end if
    end if
    if (mod(isweep, restartperiod)==0) then
      call makerestartpoint()
    end if
    call runproductiontasks
  end do
end subroutine 


!> Subroutine that gathers all the tasks related to making a restart point. 
!! These include opening (and closing) the appropriate files for molecule
!! configurations, simulation parameters and the random number generator 
!! state.
!!
!! @pre simulation state has to be initialized with mce_init.
!! @post restart files have been updated on disk or warning messages
!! have been written to standard output if something fails.
!!
subroutine makerestartpoint
  type(json_value), pointer :: parameters_json
  real(dp), allocatable :: ratios(:, :)
  parameters_json => null()
  call mce_to_json(parameters_json)
  if (replica_id == 0 .and. n_replicas > 1) then
     ! Write beta_exchange statistics.
     call be_get_acceptance_ratios(ratios)
     call json_add(parameters_json, 'beta_exchange_ratios', reshape(ratios, &
          [size(ratios)]))
     call be_reset_counters
  end if
  call json_print(parameters_json, fn_parameters_restart)
  call json_destroy(parameters_json)
end subroutine




!> All actions done during both the equilibration and production sweeps
!! should be gathered inside this routine for clarity. 
subroutine runproductiontasks
  integer :: ios
  type(json_value), pointer :: output_json
  integer :: u_output_json
  if (mod(isweep, productionperiod) == 0) then
    
    !! Write json output
    u_output_json = fileunit_getfreeunit()
    open(UNIT=u_output_json, file=fn_parameters_out, action='WRITE', &
         position='APPEND', DELIM='QUOTE', IOSTAT=ios)
    if (ios /= 0) then
       write(error_unit, *) 'ERROR: runproductiontasks failed opening ', &
            fn_parameters_out
       stop 'Program stopped by runproductiontasks.'
    end if
    call mce_to_json(output_json)
    call json_print(output_json, u_output_json)
    close(u_output_json)
    call json_destroy(output_json)
    call resetcounters
  end if
end subroutine

!> Resets the counters that are used to monitor acceptances of trial
!! moves. 
subroutine resetcounters
  call nvt_engine_reset_counters()
  if (ensemble == 'npt') then
     call npt_engine_reset_counters()
  end if
end subroutine resetcounters

!> Adjusts the maximum values for trial moves of particles and trial
!! scalings of the simulation volume. Should be used only during
!! equilibration run. 
subroutine updatemaxvalues
  !! Adjust scaling
  if (ensemble == 'npt') then
     call npt_engine_update_max_scaling()
  end if
  call nvt_engine_update_max_moves()
end subroutine updatemaxvalues
  
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
end subroutine


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

end module mc_engine
