!> Implements the streering of the simulation and controlling it's
!! input and output.
module mc_engine
use nrtype
use utils
use class_factory
use mt_stream
use particle, only: particledat
use class_poly_box
use mc_sweep
use m_fileunit
use class_parameterizer
use class_parameter_writer
use beta_exchange, only: write_stats, reset_counters
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

!> If adjusting of temperatures in the temperature series for parallel 
!! tempering is used defines the period for these actions. (DEPRECATED)
integer, save :: ptadjustperiod = 100

!> The period of writing a restart point on the disk. This includes the
!! configuration of molecules, the simulation parameters and the random
!! number generator state. 
integer, save :: restartperiod = 10000

!> The sweep counter.
integer, save :: isweep = 0

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
  
  type(particledat), dimension(:), allocatable :: particles
  !! is the pointer to the array where the particles are stored
  !! throughout the simulation.

  type(poly_box) :: simbox
  !! stores the simulation box that is used. 

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
  call getparameter(parameterreader, 'pt_adjusting_period', ptadjustperiod)
  call getparameter(parameterreader, 'restartperiod', restartperiod)

  !! Read geometry
  coordinateunit = fileunit_getfreeunit()
  statefile = 'inputconfiguration.'//trim(adjustl(idchar))
  open(file=statefile, unit=coordinateunit, action='READ', status='OLD',&
       iostat=ios)
  call factory_readstate(coordinatereader, coordinateunit, simbox, particles,&
       ios)
  if (0 /= ios) then 
    write(*, *) 'Error ', ios,' reading ', statefile, ' Stopping.' 
    stop
  end if
  close(coordinateunit)

  !! Initialize modules. 
  call mcsweep_init(parameterreader, simbox, particles)
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
  do i = 0, size(mts) - 1
     call delete(mts(i))
  end do
  if (allocated(mts)) deallocate(mts)
  if (id == 0) write (*, *) 'Program ptgbcyl was finalized succesfully.'
end subroutine 

  
!> Writes the parameters and observables of this module and its children
!! with the class_parameter_writer module and by calling the write routines
!! of the child modules.
!!
!! @param writer is the object responsible for writing the parameters.
!!
subroutine mce_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call writecomment(writer, 'mc engine parameters')
  call writeparameter(writer, 'n_equilibration_sweeps', &
  nequilibrationsweeps)
  call writeparameter(writer, 'n_production_sweeps', nproductionsweeps)
  call writeparameter(writer, 'i_sweep', isweep)
  call writeparameter(writer, 'production_period', productionperiod)
  call writeparameter(writer, 'move_adjusting_period', moveadjustperiod)
  call writeparameter(writer, 'pt_adjusting_period', ptadjustperiod)
  call writeparameter(writer, 'restartperiod', restartperiod)
  call writeparameter(writer, 'seed', seed)
  call mc_sweep_writeparameters(writer)
end subroutine


!> Runs the simulation. 
subroutine run
  do while (isweep < nequilibrationsweeps + nproductionsweeps)
    isweep = isweep + 1
    call sweep(mts, isweep)
    if (isweep <= nequilibrationsweeps) then
      call runequilibrationtasks
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
  type(particledat), allocatable :: particles(:)
  type(poly_box) :: simbox

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

  call get_system(simbox, particles)
  call factory_writestate(restartwriter, configurationunit, simbox, particles)
  close(configurationunit)
  if (allocated(particles)) deallocate(particles)

end subroutine


!> All actions performed only during equilibration sweeps and not
!! during production sweeps are gathered inside this routine.
subroutine runequilibrationtasks
  if (moveadjustperiod /= 0) then
    if (mod(isweep, moveadjustperiod) == 0) then
      call updatemaxvalues
    end if
  end if
end subroutine 

!> All actions done during both the equilibration and production sweeps
!! should be gathered inside this routine for clarity. 
subroutine runproductiontasks
  integer :: ios
  character(len=80) :: parameterfile
  type(parameter_writer) :: writer
  type(factory) :: coordinatewriter
  type(poly_box) :: simbox
  type(particledat), allocatable :: particles(:)
  integer :: be_unit
  if (mod(isweep, productionperiod) == 0) then
    !! Record snapshot of molecules and geometry.
    call get_system(simbox, particles)
    call factory_writestate(coordinatewriter, coordinateunit, simbox, particles)
    if (allocated(particles)) deallocate(particles) 
    
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
  
end module
