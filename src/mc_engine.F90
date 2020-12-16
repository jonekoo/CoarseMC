!> Implements the steering of the simulation and controlling it's
!! input and output.
module mc_engine
  !$ use omp_lib
  use m_npt_engine, only: npt_engine_reset_counters, &
       npt_engine_update_max_scaling, npt_engine_finalize, &
       npt_engine_to_json, npt_engine_init_json, pressure, npt_update
  use m_nvt_engine, only: temperature, nvt_update, &
       nvt_engine_update_max_moves, nvt_engine_reset_counters, &
       nvt_engine_to_json, nvt_engine_init_json, total_energy, simbox, &
       nvt_engine_update_decomposition, check_simbox
  use mt_stream
  use m_fileunit
  use class_poly_box, only: poly_box, getx, gety, getz
  use m_point, only: point, pair_interaction_ptr, &
       particlearray_wrapper, single_interaction_ptr
  use particle_mover, only: particlemover_init,&
       get_max_translation, getmaxmoves, setmaxmoves, particlemover_to_json
  use beta_exchange, only: write_stats, be_reset_counters, &
       beta_exchange_init => init, be_finalize => finalize, try_beta_exchanges,&
       be_get_acceptance_ratios
  use num_kind, only: dp
  use iso_fortran_env, only: error_unit, output_unit
  use json_module
  use m_json_wrapper, only: get_parameter
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
    !! Initialize modules.
    call particlemover_init(json_val)
    call nvt_engine_init_json(json_val)
    if (ensemble == 'npt') then
       call npt_engine_init_json(json_val)
    end if
    
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
    call total_energy(etotal, err)
    if (err /= 0) then
       write(error_unit, *) &
            'ERROR: mce_init: total_energy returned err = ', err
       stop 'Program stopped by mce_init_common.'
    end if
    call makerestartpoint
  end subroutine mce_init_common
  

subroutine mce_from_json(json_val)
  type(json_value), pointer, intent(in) :: json_val
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
  ensemble = 'npt'
  call get_parameter(json_val, 'ensemble', ensemble)
end subroutine mce_from_json

!> Initializes independent random number generator states for each MPI
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



!> Finalizes the simulation.
!!
!! @param id is the MPI process id.
!!
subroutine finalize(id)
  integer, intent(in) :: id
  integer :: i
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
  if (allocated(mts)) deallocate(mts)
  if (id == 0) write (output_unit, *) &
       'Program coarsemc was finalized succesfully.'
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
subroutine sweep(genstates, isweep)
  type(mt_state), intent(inout) :: genstates(0:)
  integer, intent(in) :: isweep
  real(dp) :: beta, dE
#ifdef DEBUG
  integer :: err
  write(output_unit, *) 'Move particles'
#endif
  dE = 0.
  call nvt_update(genstates, dE)
#ifdef DEBUG
  write(output_unit, *) "etotal + dE = ", etotal, " + ", dE, " = ", etotal + dE
#endif
  etotal = etotal + dE
#ifdef DEBUG
  call total_energy(etotal, err)
  write(output_unit, *) "etotal updated = ", etotal
#endif
  if (ensemble == 'npt') then
#ifdef DEBUG
     write(output_unit, *) 'Make volume move(s)'
#endif
     call npt_update(genstates(0), etotal)
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
  call json_create_object(json_val, 'mc_engine')
  call json_add(json_val, 'replica_id', replica_id)
  call json_add(json_val, 'n_equilibration_sweeps', nequilibrationsweeps)
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

  call particlemover_to_json(json_val)
  call nvt_engine_to_json(json_val)
  if (ensemble == 'npt') then
     call npt_engine_to_json(json_val)
  end if

end subroutine mce_to_json


!> Runs the simulation. 
subroutine run
  integer :: i
  do while (isweep < nequilibrationsweeps + nproductionsweeps)
    isweep = isweep + 1
    call sweep(mts, isweep)
    if (isweep <= nequilibrationsweeps) then
      if (moveadjustperiod /= 0) then
         if (mod(isweep, moveadjustperiod) == 0) then
            call updatemaxvalues
            call nvt_engine_update_decomposition
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
  
end module mc_engine
