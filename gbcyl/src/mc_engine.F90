!> Implements the streering of the simulation and controlling it's
!! input and output.
module mc_engine
  use m_particlegroup, only: make_particle_moves, update_volume, &
       total_energy, mc_sweep_writeparameters, mcsweep_finalize, &
       get_maxscaling, set_maxscaling, particlegroup, particlegroup_ptr, &
       mc_sweep_to_json, particlegroup_init
  use m_particle_factory, only: factory, factory_readstate, &
       factory_writestate, read_group
  use mt_stream
  use m_fileunit
  use class_poly_box, only: poly_box, volume, getx, gety, getz
  use class_parameterizer
  use class_parameter_writer
  use particle, only: pair_interaction, particledat, pair_interaction_ptr, &
       particlearray_wrapper
  use particle_mover, only: particlemover_init, particlemover_writeparameters,&
       get_max_translation, getmaxmoves, setmaxmoves, particlemover_to_json
  use beta_exchange, only: write_stats, reset_counters, &
       beta_exchange_init => init, be_finalize => finalize, try_beta_exchanges
  use class_pair_potential, only: conditional_pair_interaction
  use num_kind
  use iso_fortran_env, only: error_unit, output_unit
  use particlewall, only: particlewall_potential, particlewall_init, &
       particlewall_writeparameters, particlewall_to_json
  use utils, only: splitstr, join
  !$ use omp_lib
  use json_module
  use m_json_wrapper, only: get_parameter
  use m_interaction_factory, only: create_pair_interaction
  use m_particlejson_parser, only: particlearray_from_json
  implicit none
  
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
  
  !> The simulation temperature.
  real(dp), save :: temperature = -1._dp

  !> The total energy.
  real(dp), save :: etotal = 0._dp
  
  !> The simulation pressure. Only meaningful in a constant-pressure
  !! simulation.
  real(dp), save :: pressure = -1._dp
  
  !> Output unit and filename for writing parameters.
  integer, save :: pwunit

  !> The input and output unit used for reading and writing the geometry of 
  !! molecules and the simulation box.
  integer, save :: coordinateunit
  
  !> The (MPI) id of this process formatted to a character.
  character(len = 9), save :: idchar

  !> Output file name for parameters.
  character(len=:), allocatable, save :: fn_parameters_out
  character(len=:), allocatable, save :: fn_parameters_restart
  character(len=:), allocatable, save :: fn_coordinates_restart
  
  !> The random number generator states.
  type(mt_state), allocatable, save :: mts(:)
  
  !> The random number generator seed.
  integer, save :: seed = -1
  
  !> Counter for trial particle moves.
  integer, save :: nmovetrials = 0
  
  !> Counter for trial volume updates.
  integer, save :: nscalingtrials = 0
  
  !> The number of accepted particle moves
  integer, save :: nacceptedmoves = 0
  
  !> The number of accepted volume moves
  integer, save :: nacceptedscalings = 0
  
  !> The desired acceptance ratio for trial particle moves.
  real(dp), save :: moveratio
  
  !> The desired acceptance ratio for trial volume updates.
  real(dp), save :: scalingratio

  type(particlegroup_ptr), allocatable, save :: groups(:)
  type(poly_box), save :: simbox

  character(len=20), allocatable, save :: group_names(:)
  type(pair_interaction_ptr), allocatable, save :: pair_interactions(:, :)

contains
  
  subroutine mce_init_json(id, n_tasks, parameter_infile, parameter_outfile, &
       parameter_restartfile)
    integer, intent(in) :: id, n_tasks    
    character(len=*), intent(in) :: parameter_infile
    character(len=*), intent(in) :: parameter_outfile
    character(len=*), intent(in) :: parameter_restartfile
    !character(len=*), intent(in) :: coordinate_restartfile
    type(json_file) :: json
    type(json_value), pointer :: json_val, box_json
    logical :: status_ok
    character(len=:), allocatable :: error_msg
    fn_parameters_out = parameter_outfile
    fn_parameters_restart = parameter_restartfile
    !fn_coordinates_restart = coordinate_restartfile
    write(idchar, '(I9)') id 
    call json_parse(parameter_infile, json_val)
    if (json_failed()) then
       call json_print_error_message(error_unit)
       stop
    end if
    !call json_check_for_errors(status_ok, error_msg)
    !if (.not. status_ok) then
    !   write(error_unit, *) 'ERROR: ' // error_msg
    !   call json_clear_exceptions()
    !end if
    !call json%get('$', json_val)
    call mce_from_json(json_val)
    if (json_failed()) then
        call json_check_for_errors(status_ok, error_msg)
        write(error_unit, *) 'Error: '//error_msg
        call json_clear_exceptions()
        stop 'json_failed'
    !    !call json%destroy()
    end if
    call create_interactions_json(json_val, group_names, pair_interactions)
    !if (json_failed()) then
    !    call json_check_for_errors(status_ok, error_msg)
    !    write(error_unit, *) 'Error: '//error_msg
    !    call json_clear_exceptions()
    !    stop 'json_failed'
    !    !call json%destroy()
    !end if
    !call json%destroy() ! Causes segmentation fault.
    !write(*, *) nequilibrationsweeps
    call json_get(json_val, 'box', box_json)
    call simbox%from_json(box_json)

    call create_groups(json_val, simbox, group_names, pair_interactions, groups)
    
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
    character(len = 50) :: statefile 
    integer :: ios
    type(factory) :: coordinatereader
    type(particledat), allocatable :: particles(:)
    integer :: err
    type(particlearray_wrapper), allocatable :: temp_groups(:)
    integer :: i
    
    call mce_init_rng(id, n_tasks)
    
    if (temperature < 0._dp) then
       write(error_unit, *) 'ERROR: mc_engine '//&
            'tried to set a negative temperature.'
       stop 'Program stopped by mce_init_common.'  
    end if
    call beta_exchange_init(1._dp / temperature)
    
    if (pressure < 0._dp) then
       write(error_unit, *) 'ERROR: mc_engine '//&
            'tried to set a negative pressure.'
       stop 'Program stopped by mce_init_common.'  
    end if

    !! statefile = 'inputconfiguration.'//trim(adjustl(idchar))

    !! !! Read geometry
    !! if (.false.) then
    !!    allocate(temp_groups(size(group_names)))
    !!    do i = 1, size(group_names)
    !!       call read_group(statefile, group_names(i), temp_groups(i))
    !!    end do
    !! else
    !!    coordinateunit = fileunit_getfreeunit()
    !!    open(file=statefile, unit=coordinateunit, action='READ', status='OLD',&
    !!         iostat=ios)
    !!    call factory_readstate(coordinatereader, coordinateunit, simbox, &
    !!         particles, ios)
    !!    if (0 /= ios) then 
    !!       write(error_unit, *) 'ERROR ', ios,' reading ', statefile
    !!       stop
    !!    end if
    !!    close(coordinateunit)
    !! end if
    !! call create_groups(simbox, particles, group_names, &
    !!     pair_interactions, groups)


    
    !! Open output for geometries
    coordinateunit = fileunit_getfreeunit()
    statefile = 'configurations.' // trim(adjustl(idchar))
    open(file=statefile, unit=coordinateunit, action='WRITE', &
         position='APPEND', status='UNKNOWN', form='formatted', iostat=ios)
    if (0 /= ios) then
       write(error_unit, *) 'ERROR: mce_init: Failed opening ', statefile, & 
            ' for writing.'
       stop 'Program stopped by mce_init_common.'
    end if
    call total_energy(groups, simbox, pair_interactions, &
         particlewall_potential, etotal, err)
    if (err /= 0) then
       write(error_unit, *) 'ERROR: mce_init: total_energy returned err = ', err
       stop 'Program stopped by mce_init_common.'
    end if
    call makerestartpoint
  end subroutine mce_init_common
  

subroutine mce_from_json(json_val)
  type(json_value), pointer, intent(in) :: json_val
  character(kind=CK, len=:), allocatable :: msg
  logical :: status_ok
  call get_parameter(json_val, 'n_equilibration_sweeps', &
       nequilibrationsweeps, error_lb=0)
  call get_parameter(json_val, 'n_production_sweeps', &
       nproductionsweeps, error_lb=0)
  call get_parameter(json_val, 'production_period', &
       productionperiod, error_lb=1)
  call get_parameter(json_val, 'i_sweep', isweep, error_lb=0, &
       warn_ub=nequilibrationsweeps + nproductionsweeps)  
  call get_parameter(json_val, 'move_adjusting_period', &
       moveadjustperiod, error_lb=1)
  call get_parameter(json_val, 'move_ratio', moveratio, error_lb=0._dp, &
       error_ub=1._dp)
  call get_parameter(json_val, 'scaling_ratio', scalingratio, error_lb=0._dp, &
       error_ub=1._dp)
  call get_parameter(json_val, 'pt_period', pt_period, error_lb=1, &
       warn_ub=nequilibrationsweeps + nproductionsweeps)
  call get_parameter(json_val, 'restartperiod', restartperiod, error_lb=1)
  call get_parameter(json_val, 'temperature', temperature, error_lb=0._dp)
  call get_parameter(json_val, 'pressure', pressure, error_lb=0._dp)
  call get_parameter(json_val, 'seed', seed, warn_lb=1000)
  call get_parameter(json_val, 'nmovetrials', nmovetrials, error_lb=0)
  call get_parameter(json_val, 'nacceptedmoves', nacceptedmoves, &
       error_ub=nmovetrials, error_lb=0)
  call get_parameter(json_val, 'nscalingtrials', nscalingtrials, error_lb=0)
  call get_parameter(json_val, 'nacceptedscalings', nacceptedscalings, &
       error_ub=nscalingtrials, error_lb=0)
  allocate(group_names(0))
  call get_parameter(json_val, 'groups', group_names)
  !! Initialize modules.
  call particlewall_init(json_val)
  call particlemover_init(json_val)
  call particlegroup_init(json_val)
end subroutine mce_from_json


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
  !call init(mts(thread_id), seed)
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


subroutine create_groups(json_val, simbox, group_names, &
     pair_interactions, groups)
  type(json_value), pointer, intent(in) :: json_val
  type(poly_box), intent(in) :: simbox
  character(len=*), intent(in) :: group_names(:)
  character(kind=CK, len=:), allocatable :: group_name
  type(pair_interaction_ptr), intent(in) :: pair_interactions(:, :)
  type(particlegroup_ptr), allocatable, intent(out) :: groups(:)
  character(len=20) :: group_type
  integer :: i
  real(dp) :: max_cutoff
  type(json_value), pointer :: groups_json, groups_json_element
  class(particledat), allocatable :: particles(:)
  max_cutoff = pair_interactions(1, 1)%ptr%get_cutoff()
  allocate(groups(size(group_names)))
  call json_get(json_val, 'particle_groups', groups_json)
  do i = 1, json_count(groups_json)
     call json_get_child(groups_json, i, groups_json_element)
     call get_parameter(groups_json_element, 'name', group_name)
     if (any(group_names == group_name)) then
        call particlearray_from_json(groups_json_element, particles)
        allocate(groups(i)%ptr, source=particlegroup(simbox, particles, &
             min_cell_length=max_cutoff + 2 * get_max_translation(), &
             min_boundary_width=2 * get_max_translation(), name=group_name))
     end if
     deallocate(particles)
  end do
end subroutine create_groups


subroutine create_groups2(simbox, particles, group_names, &
     pair_interactions, groups)
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: particles(:)
  character(len=*), intent(in) :: group_names(:)
  type(pair_interaction_ptr), intent(in) :: pair_interactions(:, :)
  type(particlegroup_ptr), allocatable, intent(out) :: groups(:)
  character(len=20) :: group_type
  integer :: i
  real(dp) :: max_cutoff
  max_cutoff = pair_interactions(1, 1)%ptr%get_cutoff()
  allocate(groups(size(group_names)))
  do i = 1, size(group_names)
     if (group_names(i) == 'gb') then
        allocate(groups(i)%ptr, source=particlegroup(simbox, &
             pack(particles, particles%rod), &
             min_cell_length=max_cutoff + &
             2 * get_max_translation(), &
             min_boundary_width=2 * get_max_translation(), name=group_names(i)))
     else if (group_names(i) == 'lj') then
        allocate(groups(i)%ptr, source=particlegroup(simbox, &
             pack(particles, .not. particles%rod), &
             min_cell_length=max_cutoff + &
             2 * get_max_translation(), &
             min_boundary_width=2 * get_max_translation(), name=group_names(i)))
     else
        write(error_unit, *) 'ERROR: unknown group type ' // &
             trim(adjustl(group_type))
        stop 'Program stopped by create_groups2.'
     end if
  end do
end subroutine create_groups2


!> Finalizes the simulation.
!!
!! @param id is the MPI process id.
!!
subroutine finalize(id)
  integer, intent(in) :: id
  integer :: i
  close(coordinateunit)
  call makerestartpoint
  call mcsweep_finalize
  call be_finalize
  do i = 0, size(mts) - 1
     call delete(mts(i))
  end do
  if (allocated(mts)) deallocate(mts)
  if (id == 0) write (*, *) 'Program ptgbcyl was finalized succesfully.'
end subroutine 

  
!> Runs one sweep of Metropolis Monte Carlo updates to the system. A
!! full Parallel tempering NPT-ensemble sweep consists of trial moves
!! of particles, trial scaling of the simulation box (barostat) and an
!! exchange of particle system coordinates with another particle system
!! (replica) in another temperature (replica exchange).
!! 
!! @param genstates random number generator states for all threads.
!! @param isweep the sweep counter.
!!  
subroutine sweep(simbox, groups, genstates, isweep)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup_ptr), intent(inout) :: groups(:)
  type(mt_state), intent(inout) :: genstates(0:)
  integer, intent(in) :: isweep
  real(dp) :: beta, dE
  integer :: n_trials, n_accepted, err
  dE = 0.
  call make_particle_moves(groups, genstates, simbox, temperature, & 
       pair_interactions, particlewall_potential, dE, n_trials, &
       n_accepted)
  write(output_unit, *) "etotal + dE = ", etotal, " + ", dE, " = ", etotal + dE
  call total_energy(groups, simbox, pair_interactions, &
       particlewall_potential, etotal, err)
  write(output_unit, *) "etotal updated = ", etotal
  nmovetrials = nmovetrials + n_trials
  nacceptedmoves = nacceptedmoves + n_accepted
  call update_volume(simbox, groups, genstates(0), pair_interactions,&
       particlewall_potential, temperature, pressure, etotal, n_trials, &
       n_accepted)
  nscalingtrials = nscalingtrials + n_trials
  nacceptedscalings = nacceptedscalings + n_accepted
  call check_simbox(simbox)
  if (mod(isweep, pt_period) == 0) then
     beta = 1._dp / temperature
     call try_beta_exchanges(beta, etotal + pressure * volume(simbox), 3, &
          genstates(0)) 
     temperature = 1._dp / beta
  end if
end subroutine sweep


subroutine mce_to_json(json_val, coordinates_json)
  type(json_value), pointer, intent(out) :: json_val
  type(json_value), pointer, intent(out), optional :: coordinates_json
  integer :: i, j
  type(json_value), pointer :: json_child, group_name, pair_ia_json, &
       pair_ia_element, group_json, group_json_element, box_json
  call json_create_object(json_val, 'mc_engine')
  
  call json_add(json_val, 'n_equilibration_sweeps', &
  nequilibrationsweeps)
  call json_add(json_val, 'n_production_sweeps', nproductionsweeps)
  call json_add(json_val, 'i_sweep', isweep)
  call json_add(json_val, 'production_period', productionperiod)
  call json_add(json_val, 'move_adjusting_period', moveadjustperiod)
  call json_add(json_val, 'pt_period', pt_period)
  call json_add(json_val, 'restartperiod', restartperiod)
  call json_add(json_val, 'seed', seed)
  call json_add(json_val, 'pressure', pressure)
  call json_add(json_val, 'temperature', temperature)
  call json_add(json_val, 'volume', volume(simbox))
  call json_add(json_val, 'total_energy', etotal)
  call json_add(json_val, 'enthalpy', etotal + &
       volume(simbox) * pressure)
  call json_add(json_val, 'move_ratio', moveratio)
  call json_add(json_val, 'scaling_ratio', scalingratio)
  call json_add(json_val, 'nmovetrials', nmovetrials)
  call json_add(json_val, 'nacceptedmoves', nacceptedmoves)
  if (nmovetrials > 0) then
     call json_add(json_val, 'current_move_ratio', &
          real(nacceptedmoves, dp)/real(nmovetrials, dp))
  else
     call json_add(json_val, 'current_move_ratio', 'nan')
  end if
  call json_add(json_val, 'nscalingtrials', nscalingtrials)
  call json_add(json_val, 'nacceptedscalings', nacceptedscalings)
  if (nscalingtrials > 0) then
     call json_add(json_val, 'current_scaling_ratio', &
          real(nacceptedscalings, dp)/real(nscalingtrials, dp))
  else
     call json_add(json_val, 'current_scaling_ratio', 'nan')
  end if
  
  !! Write group names.
  call json_create_array(json_child, 'groups')
  do i = 1, size(group_names)
     call json_create_string(group_name, group_names(i), '')
     call json_add(json_child, group_name)
  end do

  call json_add(json_val, json_child)
  
  call json_create_array(pair_ia_json, 'pair_interactions')
  do j = 1, size(pair_interactions, 2)
     do i = 1, j
        call json_create_object(pair_ia_element, '')
        call json_add(pair_ia_element, 'participants', [group_names(i), group_names(j)])
        call pair_interactions(i, j)%ptr%to_json(pair_ia_element)
        call json_add(pair_ia_json, pair_ia_element)
     end do
  end do
  call json_add(json_val, pair_ia_json) 
  call particlewall_to_json(json_val)
  call particlemover_to_json(json_val)
  call mc_sweep_to_json(json_val)

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
  if (present(coordinates_json)) then
     call json_add(coordinates_json, group_json)
  else
     call json_add(json_val, group_json)
  end if
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
  parameters_json => null()
  call mce_to_json(parameters_json) !, coordinates_json)
  call json_print(parameters_json, fn_parameters_restart)
  call json_destroy(parameters_json)
end subroutine


subroutine writestate(writer, unit, simbox, groups)
  type(factory), intent(in) :: writer
  integer, intent(in) :: unit
  type(poly_box), intent(in) :: simbox
  type(particlegroup_ptr), intent(in) :: groups(:)
  type(particledat), allocatable :: particles(:)
  integer :: n, i
  n = 0
  do i = 1, size(groups)
     n = n + size(groups(i)%ptr%particles)
  end do
  allocate(particles(n))
  n = 0
  do i = 1, size(groups)
     particles(n + 1:n + size(groups(i)%ptr%particles)) = &
          groups(i)%ptr%particles
     n = n + size(groups(i)%ptr%particles)
  end do
  call factory_writestate(writer, unit, simbox, particles)  
end subroutine writestate


!> All actions done during both the equilibration and production sweeps
!! should be gathered inside this routine for clarity. 
subroutine runproductiontasks
  integer :: ios
  integer :: be_unit
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
        
    !! Write beta_exchange statistics:
    if (trim(adjustl(idchar)) == "0") then
      be_unit = fileunit_getfreeunit()
      open(unit=be_unit, file="beta_exchange.stats", action="WRITE", &
           position="APPEND", iostat=ios)
      if (ios /= 0) then
         write(error_unit, *) 'Warning: runproductiontasks failed opening ' //&
              'beta_exchange.stats'
      else
         call write_stats(be_unit)
         close(be_unit)
      end if
      call reset_counters
    end if
    close(pwunit)
    call resetcounters
  end if
end subroutine

!> Resets the counters that are used to monitor acceptances of trial
!! moves. 
subroutine resetcounters
  nacceptedmoves = 0
  nmovetrials = 0
  nscalingtrials = 0 
  nacceptedscalings = 0
end subroutine resetcounters

!> Adjusts the maximum values for trial moves of particles and trial
!! scalings of the simulation volume. Should be used only during
!! equilibration run. 
subroutine updatemaxvalues
  real(dp) :: newdthetamax
  real(dp) :: newdximax
  !! Adjust scaling
  if (nscalingtrials > 0) then
     call set_maxscaling(newmaxvalue(nscalingtrials, nacceptedscalings, &
          scalingratio, get_maxscaling()))
  end if
  if (nmovetrials > 0) then
     call getmaxmoves(newdximax, newdthetamax)
     !! Adjust translation
     newdximax = newmaxvalue(nmovetrials, nacceptedmoves, moveratio, newdximax)
     
     !! Update the minimum cell side length of the cell list because the
     !! maximum translation has changed: 
     
     !! This should adjust rotations < pi/2 to move the particle end as much as
     !! a random translation. 4.4 is the assumed molecule length-to-breadth 
     !! ratio.
     newdthetamax = 2._dp * asin(newdximax/4.4_dp) 
     call setmaxmoves(newdximax, newdthetamax)
  end if
end subroutine updatemaxvalues
  
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


subroutine create_interactions_json(json_val, group_names, pair_ias)
  type(json_value), pointer, intent(in) :: json_val
  character(len=*), intent(in) :: group_names(:)
  type(pair_interaction_ptr), intent(inout), allocatable :: pair_ias(:, :)
  type(json_value), pointer :: pair_ia_json, pair_ia_element
  character(len=len(group_names)), allocatable :: participants(:)
  logical :: found
  integer :: i, j, k, l, astat
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
                      'Warning: only the first interaction with participants ', &
                      participants, ' is used.'
              else
                 allocate(pair_ias(k, l)%ptr, &
                      source=create_pair_interaction(pair_ia_element), &
                      stat=astat)
                 if (astat /= 0) then
                    write(error_unit, *) 'create_interactions_json could ' // &
                         'not allocate interaction between ', group_names(k), &
                         ' and ', group_names(l), '.'
                    stop 'create_interactions_json unable to continue.'
                 end if
                 if (k /= l) then
                    pair_ias(l, k)%ptr => pair_ias(k, l)%ptr
                 end if
              end if
           else
              write(error_unit, *) 'ERROR: Pair interaction has invalid participant.'
              stop 'create_interactions_json unable to continue.'
           end if
        else
           write(error_unit, *) 'ERROR: Pair interaction has ', &
                size(participants), ' participants.'
           stop 'create_interactions_json unable to continue.'
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
           write(error_unit, *) 'Interaction between ', group_names(i), " and ", &
                group_names(j), 'is not defined.'
           stop 'create_interactions_json unable to continue.'
        end if
     end do
  end do
end subroutine create_interactions_json

end module mc_engine
