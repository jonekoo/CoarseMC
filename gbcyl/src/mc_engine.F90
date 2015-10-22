!> Implements the streering of the simulation and controlling it's
!! input and output.
module mc_engine
  use m_particlegroup, only: make_particle_moves, update_volume, get_total_energy, &
       mc_sweep_writeparameters, mcsweep_finalize, mcsweep_init, &
       get_maxscaling, set_maxscaling, particlegroup
  
  use class_factory, only: factory, factory_readstate, factory_writestate
  use mt_stream
  use m_fileunit
  use class_poly_box, only: poly_box, volume, getx, gety, getz
  use class_parameterizer
  use class_parameter_writer
  use particle, only: pair_interaction
  use particle_mover, only: particlemover_init, particlemover_writeparameters,&
       get_max_translation, getmaxmoves, setmaxmoves
  use beta_exchange, only: write_stats, reset_counters, &
       beta_exchange_init => init, be_finalize => finalize, try_beta_exchanges
  use class_pair_potential, only: create_conditional_interaction, pp_init, &
       pp_writeparameters
  use iso_fortran_env, only: dp => REAL64, error_unit
  use particlewall, only: particlewall_potential, particlewall_init, &
       particlewall_writeparameters
  !$ use omp_lib
  implicit none
  private
  
  public :: mce_init
  public :: run
  public :: finalize
  public :: mce_writeparameters
  
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
  
  type(particlegroup), save :: group
  type(poly_box), save :: simbox

  class(pair_interaction), allocatable :: pair_ia
  
contains
  
!> Initializes this module and its dependencies.
!!
!! @param id the process id for an MPI process.
!! @param n_tasks total number of MPI processes.
!!
subroutine mce_init(id, n_tasks)
  integer, intent(in) :: id
  integer, intent(in) :: n_tasks 
  logical :: isrestart 
  character(len = 50) :: statefile 
  integer :: ios
  type(factory) :: coordinatereader
  character(len=80) :: parameterinputfile

  type(parameterizer) :: parameterreader
  integer :: thread_id = 0, n_threads = 1
  
  write(idchar, '(I9)') id 
  idchar = trim(adjustl(idchar))
  isrestart = .false.
  parameterinputfile='inputparameters.'//trim(adjustl(idchar))
  parameterreader = new_parameterizer(parameterinputfile, iostat=ios)
  if (ios/=0) then
    write(*,*) 'Could not open ' // trim(adjustl(parameterinputfile)) // &
         '. Stopping.'
    stop
  end if
  !$ n_threads = omp_get_max_threads()
  allocate(mts(0:n_threads - 1))
  call set_mt19937
  call new(mts(thread_id))
  call getparameter(parameterreader, 'seed', seed)
  if (seed < 0) then 
    call system_clock(seed)
    write(*, *) "Seeding RNG with system clock."
  end if
  if (seed < 0) then
    write(*, *) "system_clock query failed, using default seed 1234567."
    seed = 1234567 
  end if
  call init(mts(thread_id), seed)

  !$ write(*, *) 'Running with ', n_threads, ' threads.'

  !! Give different random number streams to each OpenMP thread inside a
  !! MPI task.
  !$OMP PARALLEL DO
  !$ do thread_id = 1, n_threads-1
    if (id + n_tasks * thread_id > 0) then 
      call create_stream(mts(0), mts(thread_id), id + n_tasks * thread_id)
    end if
  !$ end do
  !$OMP END PARALLEL DO

  call getparameter(parameterreader, 'n_equilibration_sweeps', &
  nequilibrationsweeps)
  call getparameter(parameterreader, 'n_production_sweeps', &
  nproductionsweeps)
  call getparameter(parameterreader, 'production_period', &
  productionperiod)
  call getparameter(parameterreader, 'i_sweep', isweep)
  call getparameter(parameterreader, 'move_adjusting_period', &
       moveadjustperiod)
  call getparameter(parameterreader, 'move_ratio', moveratio)
  call getparameter(parameterreader, 'scaling_ratio', scalingratio)
  call getparameter(parameterreader, 'pt_period', pt_period)
  call getparameter(parameterreader, 'pt_adjusting_period', ptadjustperiod)
  call getparameter(parameterreader, 'restartperiod', restartperiod)
  call getparameter(parameterreader, 'temperature', temperature)
  if (temperature < 0._dp) then
     write(error_unit, *) 'mc_engine:'//&
          'trying to set a negative temperature, stopping.'
     stop  
  end if
  call getparameter(parameterreader, 'pressure', pressure)
  if (pressure < 0._dp) then
     write(error_unit, *) 'mc_engine:'//&
          'trying to set a negative pressure, stopping.'
     stop  
  end if

  call getparameter(parameterreader, 'nmovetrials', nmovetrials)
  call getparameter(parameterreader, 'nscalingtrials', nscalingtrials)
  call getparameter(parameterreader, 'nacceptedscalings', nacceptedscalings)

  !! Read geometry
  coordinateunit = fileunit_getfreeunit()
  statefile = 'inputconfiguration.'//trim(adjustl(idchar))
  open(file=statefile, unit=coordinateunit, action='READ', status='OLD',&
       iostat=ios)
  call factory_readstate(coordinatereader, coordinateunit, simbox, &
       group%particles, ios)
  if (0 /= ios) then 
    write(*, *) 'Error ', ios,' reading ', statefile, ' Stopping.' 
    stop
  end if
  close(coordinateunit)

  !! Initialize modules.
  !call energy_init(parameterreader)
  call pp_init(parameterreader)
  allocate(pair_ia, source=create_conditional_interaction())
  call particlewall_init(parameterreader)
  call particlemover_init(parameterreader)
  call mcsweep_init(parameterreader)
  call beta_exchange_init(1._dp / temperature)
  call group%init(simbox, min_cell_length=pair_ia%get_cutoff() + &
       2 * get_max_translation(), min_boundary_width=2 * get_max_translation())
  call delete(parameterreader)
 
  !! Open output for geometries
  coordinateunit = fileunit_getfreeunit()
  statefile = 'configurations.' // trim(adjustl(idchar))
  open(file=statefile, unit=coordinateunit, action='WRITE', &
       position='APPEND', status='UNKNOWN', form='formatted', iostat=ios)
  if (0 /= ios) then
    write(*, *) 'mce_init: Failed opening ', statefile, & 
         ' for writing. Stopping.'
    stop
  end if
  call makerestartpoint
end subroutine 


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
subroutine sweep(simbox, group, genstates, isweep)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup), intent(inout) :: group
  type(mt_state), intent(inout) :: genstates(0:)
  integer, intent(in) :: isweep
  real(dp) :: beta
  integer :: n_trials, n_accepted
  call make_particle_moves(simbox, group, genstates, pair_ia, &
       particlewall_potential, temperature, n_trials, n_accepted)
  nmovetrials = nmovetrials + n_trials
  nacceptedmoves = nacceptedmoves + n_accepted

  call update_volume(simbox, group, genstates(0), pair_ia, &
       particlewall_potential, temperature, pressure, &
       n_trials, n_accepted)
  nscalingtrials = nscalingtrials + n_trials
  nacceptedscalings = nacceptedscalings + n_accepted
  call check_simbox(simbox)
  if (mod(isweep, pt_period) == 0) then
     beta = 1._dp / temperature
     call try_beta_exchanges(beta, get_total_energy(), 3, genstates(0)) 
     temperature = 1._dp / beta
  end if
end subroutine sweep

!> Writes the parameters and observables of this module and its children
!! with the class_parameter_writer module and by calling the write routines
!! of the child modules.
!!
!! @param writer is the object responsible for writing the parameters.
!!
subroutine mce_writeparameters(writer)
  type(parameter_writer), intent(inout) :: writer
  call writecomment(writer, 'mc engine parameters')
  call writeparameter(writer, 'n_equilibration_sweeps', &
  nequilibrationsweeps)
  call writeparameter(writer, 'n_production_sweeps', nproductionsweeps)
  call writeparameter(writer, 'i_sweep', isweep)
  call writeparameter(writer, 'production_period', productionperiod)
  call writeparameter(writer, 'move_adjusting_period', moveadjustperiod)
  call writeparameter(writer, 'pt_period', pt_period)
  call writeparameter(writer, 'pt_adjusting_period', ptadjustperiod)
  call writeparameter(writer, 'restartperiod', restartperiod)
  call writeparameter(writer, 'seed', seed)
  call writeparameter(writer, 'pressure', pressure)
  call writeparameter(writer, 'temperature', temperature)
  call writeparameter(writer, 'volume', volume(simbox))
  call writeparameter(writer, 'enthalpy', get_total_energy() + &
       volume(simbox) * pressure)
  call writeparameter(writer, 'move_ratio', moveratio)
  call writeparameter(writer, 'scaling_ratio', scalingratio)
  call writeparameter(writer, 'nmovetrials', nmovetrials)
  call writeparameter(writer, 'nacceptedmoves', nacceptedmoves)
  if (nmovetrials > 0) then
     call writeparameter(writer, 'current_move_ratio', &
          real(nacceptedmoves, dp)/real(nmovetrials, dp))
  else
     call writeparameter(writer, 'current_move_ratio', 'nan')
  end if
  call writeparameter(writer, 'nscalingtrials', nscalingtrials)
  call writeparameter(writer, 'nacceptedscalings', nacceptedscalings)
  if (nscalingtrials > 0) then
     call writeparameter(writer, 'current_scaling_ratio', &
          real(nacceptedscalings, dp)/real(nscalingtrials, dp))
  else
     call writeparameter(writer, 'current_scaling_ratio', 'nan')
  end if
  !call energy_writeparameters(writer)
  call pp_writeparameters(writer)
  call particlewall_writeparameters(writer)
  call particlemover_writeparameters(writer)
  call mc_sweep_writeparameters(writer)
end subroutine


!> Runs the simulation. 
subroutine run
  do while (isweep < nequilibrationsweeps + nproductionsweeps)
    isweep = isweep + 1
    call sweep(simbox, group, mts, isweep)
    if (isweep <= nequilibrationsweeps) then
      if (moveadjustperiod /= 0) then
         if (mod(isweep, moveadjustperiod) == 0) then
            call updatemaxvalues
            group%sl%min_length = pair_ia%get_cutoff() + 2._dp * get_max_translation()
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
  type(factory) :: restartwriter
  integer :: configurationunit
  integer :: parameterunit
  type(parameter_writer) :: pwriter
  character(len=80) :: parameterfile
  character(len=80) :: configurationfile
  integer :: ios

  !! Write parameters to a restartfile
  parameterunit = fileunit_getfreeunit()
  parameterfile = 'restartparameters.'//idchar
  open(UNIT=parameterunit, FILE=parameterfile, action='WRITE',&
       status='REPLACE', DELIM='QUOTE', iostat=ios)
  if (ios /= 0) write(*, *) 'makerestartpoint: Warning: Failed opening', &
  parameterfile
  pwriter = new_parameter_writer(parameterunit)
  call mce_writeparameters(pwriter)
  close(parameterunit)

  !! Write configurations to a restartfile
  configurationunit = fileunit_getfreeunit()
  configurationfile = 'restartconfiguration.'//idchar
  open(file=configurationfile, unit=configurationunit, &
       action='WRITE', status='REPLACE', iostat=ios)
  if (ios /= 0) write(*, *) 'makerestartpoint: Warning: Failed opening', &
       configurationfile

  call factory_writestate(restartwriter, configurationunit, simbox, &
       group%particles)
  close(configurationunit)

end subroutine

!> All actions done during both the equilibration and production sweeps
!! should be gathered inside this routine for clarity. 
subroutine runproductiontasks
  integer :: ios
  character(len=80) :: parameterfile
  type(parameter_writer) :: writer
  type(factory) :: coordinatewriter
  integer :: be_unit
  if (mod(isweep, productionperiod) == 0) then
    !! Record snapshot of molecules and geometry.
    call factory_writestate(coordinatewriter, coordinateunit, simbox, &
         group%particles)
    
    !! Record simulation parameters.
    parameterfile = 'parameters.'//trim(adjustl(idchar))
    pwunit = fileunit_getfreeunit()
    open(UNIT=pwunit, FILE=parameterfile, action='WRITE', position='APPEND',&
         DELIM='QUOTE', IOSTAT=ios) 
    if (ios /= 0) then
      write(*, *) 'runproductiontasks: error opening', parameterfile
      stop
    end if
    writer = new_parameter_writer(pwunit)
    call mce_writeparameters(writer)

    !! Write beta_exchange statistics:
    if (trim(adjustl(idchar)) == "0") then
      be_unit = fileunit_getfreeunit()
      open(unit=be_unit, file="beta_exchange.stats", action="WRITE", &
           position="APPEND")
      call write_stats(be_unit)
      close(be_unit)
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
  call set_maxscaling(newmaxvalue(nscalingtrials, nacceptedscalings, &
       scalingratio, get_maxscaling()))
  call getmaxmoves(newdximax, newdthetamax)
  !! Adjust translation
  newdximax = newmaxvalue(nmovetrials, nacceptedmoves, moveratio, newdximax)
  
  !! Update the minimum cell side length of the cell list because the maximum
  !! translation has changed: 
  
  !! This should adjust rotations < pi/2 to move the particle end as much as
  !! a random translation. 4.4 is the assumed molecule length-to-breadth 
  !! ratio.
  newdthetamax = 2._dp * asin(newdximax/4.4_dp) 
  call setmaxmoves(newdximax, newdthetamax)
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
  if (simbox%xperiodic .and. getx(simbox) < 2._dp * pair_ia%get_cutoff()) &
       stop 'Simulation box too small!'
  if (simbox%yperiodic .and. gety(simbox) < 2._dp * pair_ia%get_cutoff()) &
       stop 'Simulation box too small!'
  if (simbox%zperiodic .and. getz(simbox) < 2._dp * pair_ia%get_cutoff()) &
       stop 'Simulation box too small!'
end subroutine

end module
