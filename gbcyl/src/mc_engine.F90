!> A module for streering the simulation and controlling it's input and output.
module mc_engine
use nrtype
use utils
use class_factory
use mt_stream
use particle, only: particledat, initParticle, particle_writeparameters
use class_poly_box
use mc_sweep, only: mc_sweep_init => init, updatemaxvalues, sweep, &
mc_sweep_writeparameters, resetcounters, gettemperature, settemperature, get_system, set_system 
use pt
use m_fileunit
use class_parameterizer
use class_parameter_writer
!$ use omp_lib
implicit none
private

public :: init
public :: run
public :: finalize
public :: mce_writeparameters

interface init
  module procedure mce_init
end interface
  
integer, save :: nequilibrationsweeps = 0
integer, save :: nproductionsweeps = 0
!! define the number of equilibration and production MC sweeps in the 
!! simulation. total no. sweeps is nequilibrationsweeps + nproductionsweeps.
!! Equilibration sweeps may not obey detailed balance. 

integer, save :: productionperiod = 1
!! defines the period of writing the configuration of molecules and simulation
!! parameters on the disk. If productionperiod=400, every 400th MC sweep the 
!! configuration and parameters are written to the disk. 

integer, save :: moveadjustperiod = 100
!! defines the period of adjusting maximum trial move sizes during the
!! equilibration MC sweeps.

integer, save :: ptadjustperiod = 100
!! if adjusting of temperatures in the temperature series for parallel 
!! tempering is used defines the period for these actions.

integer, save :: restartperiod = 10000
!! defines the period of writing a restart point on the disk. This includes the
!! configuration of molecules, the simulation parameters and the random number
!! generator state. 

integer, save :: isweep = 0
!! the sweep counter.

integer, save :: pwunit
character(len=80), save :: parameterfilename
!! output unit and filename for the para

integer, save :: coordinateunit
!! The input and output unit used for reading and writing the geometry of 
!! molecules and the simulation box.

character(len = 9), save :: idchar
!! The (MPI) id of this process formatted to a character.

!logical :: isopenmp = .false.
type(mt_state), allocatable, save :: mts(:)
type(mt_state), save :: mts_new

integer, save :: seed = 123456

contains
  
!> Initializes the state of the simulation. 
!!
!! pre-condition: state must not be initialized by other means. 
!! pre-condition: MPI_INIT has been called.
!! post-condition: simulation can be run.  
!!
!! @p the process id for an MPI process.
!! @p ispt True if a parallel tempering simulation is initialized.
!!
subroutine mce_init(id, n_tasks)
  integer, intent(in) :: id
  integer, intent(in) :: n_tasks 
  logical :: isrestart 
  character(len = 50) :: statefile 
  integer :: ios
  type(factory) :: coordinatereader
  character(len=80) :: parameterinputfile

  integer, save :: rngunit
  character(len = 80), save :: rngfile = 'mtstate.'
  !! output unit and file for saving the random number generator state for a 
  !! restart.

  type(parameterizer) :: parameterreader
  integer :: thread_id = 0, n_threads = 1
  
  type(particledat), dimension(:), allocatable :: particles
  !! is the pointer to the array where the particles are stored throughout the
  !! simulation.

  type(poly_box) :: simbox
  !! stores the simulation box that is used. 

  write(idchar, '(I9)') id 
  idchar = trim(adjustl(idchar))
  isrestart = .false.
  parameterinputfile='inputparameters.'//trim(adjustl(idchar))
  parameterreader = new_parameterizer(parameterinputfile, iostat=ios)
  if (ios/=0) then
    write(*,*) 'Could not open ', parameterinputfile, '. Stopping.'
    stop
  end if
  rngfile = trim(adjustl(rngfile)) // trim(adjustl(idchar))
  rngunit = fileunit_getfreeunit()
  !! For mtmod use
  !!open(unit = rngunit, file = trim(adjustl(rngfile)), &
  !!action = 'READ', status = 'OLD', form = 'FORMATTED', iostat = ios)
  !! For mt_stream use
  open(unit = rngunit, file = trim(adjustl(rngfile)), &
  action = 'READ', status = 'OLD', form = 'UNFORMATTED', iostat = ios)
  allocate(mts(0:0))
  call set_mt19937
  call new(mts(thread_id))
  if(0 /= ios) then
    write(*, *) 'Warning: Could not open ' // trim(rngfile) // &
    ' for reading!'
    call getparameter(parameterreader, 'seed', seed)
    !call sgrnd(seed + id)
    call init(mts(thread_id), seed)
  else
    !call mtget(rngunit, 'f')
    call read(mts(thread_id), rngunit)
    close(rngunit)    
  end if

  !$ n_threads = omp_get_max_threads()
  !$ write(*, *) 'Running with ', n_threads, ' threads.'
  mts_new = mts(0)
  !$ if (allocated(mts)) deallocate(mts)
  !$ allocate(mts(0:n_threads-1))
  !! Give different random number streams to each OpenMP thread inside a MPI task

  !$ mts(0) = mts_new
  !! The structure below will probably make a mess of the streams read from files.
  !$ do thread_id = 1, n_threads-1
    if (id + n_tasks * thread_id > 0) then 
      call create_stream(mts_new, mts(thread_id), id + n_tasks * thread_id)
    end if
  !$ end do

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
!  call getparameter(parameterreader, 'isopenmp', isopenmp)
  call getparameter(parameterreader, 'restartperiod', restartperiod)

  !! Read geometry
  coordinateunit = fileunit_getfreeunit()
  statefile = 'inputconfiguration.'//trim(adjustl(idchar))
  !! In any case the output should be appended to configurations.(id) and restartfile
  !! should be used only for reading once and then overwriting the file.
  open(file=statefile, unit=coordinateunit, action='READ', status='OLD', iostat=ios)
  call readstate(coordinatereader, coordinateunit, simbox, particles, ios)
  if (0/=ios) then 
    write(*, *) 'Error ', ios,' reading ', statefile, ' Stopping.' 
    stop
  end if
  close(coordinateunit)

  !! Initialize modules. 
  call initparticle(parameterreader)
!  if (isopenmp) then 
!    call ompsweep_init(parameterreader, simbox, particles)
!  else
    call mc_sweep_init(parameterreader, simbox, particles)
!  end if
  call delete(parameterreader)
 
  !! Open output for geometries
  coordinateunit = fileunit_getfreeunit()
  statefile='configurations.'//trim(adjustl(idchar))
  open(file=statefile,unit=coordinateunit,action='WRITE',position='APPEND',&
  status='UNKNOWN',form='formatted',iostat=ios)
  if (0/=ios) then
    write(*,*) 'mce_init: Failed opening ',statefile, ' for writing. Stopping.'
    stop
  end if
  call makerestartpoint
end subroutine 

!> Writes the parameters and observables of this module and its children
!! with the class_parameter_writer module and by calling the write routines
!! of the child modules.
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
!  call writeparameter(writer, 'isopenmp', isopenmp)
  call writeparameter(writer, 'restartperiod', restartperiod)
  call writeparameter(writer, 'seed', seed)
!  if (isopenmp) then
!    call ompsweep_writeparameters(writer)
!  else
    call mc_sweep_writeparameters(writer)
!  end if
  call particle_writeparameters(writer)
end subroutine

!> Finalizes the simulation.
!!
!! pre-condition: simulation must be initialized.
!! post-conditions: 1. memory allocated by this program is freed.
!!
subroutine finalize(id)
  integer, intent(in) :: id
  close(coordinateunit)
  call makerestartpoint
  !if (associated(particles)) deallocate(particles)
  if (id == 0) write (*, *) 'End of program gbcyl.'
  call pt_finalize()
end subroutine 
  
!> Runs the simulation. 
!! 
!! pre-condition: simulation state has to be initialized. 
!! 
subroutine run
  do while (isweep < nequilibrationsweeps + nproductionsweeps)
    isweep = isweep + 1
!    if (isopenmp) then 
!      call ompsweep(simbox, particles, mts)
!    else
      call sweep(mts, isweep)
!    end if
    if (isweep <= nequilibrationsweeps) then
      call runequilibrationtasks
    end if
    call runproductiontasks
    if (mod(isweep, restartperiod)==0) then
      call makerestartpoint()
    end if
  end do
end subroutine 


!> Subroutine that gathers all the tasks related to making a restart point. 
!! These include opening (and closing) the appropriate files for molecule
!! configurations, simulation parameters and the random number generator 
!! state.
!!
!! pre-condition: simulation state has to be initialized with mce_init.
!! post-condition: restart files have been updated on disk or warning messages
!! have been written to standard output if something fails.
!!
subroutine makerestartpoint
  type(factory) :: restartwriter
  integer :: configurationunit
  integer :: parameterunit
  type(parameter_writer) :: pwriter
  character(len=80) :: parameterfile
  character(len=80) :: configurationfile

  integer :: rngunit
  character(len=80) :: rngfile

  integer :: ios
  type(particledat), allocatable :: particles(:)
  type(poly_box) :: simbox


  !! Write parameters to a restartfile
  parameterunit = fileunit_getfreeunit()
  open(UNIT=parameterunit,FILE='restartparameters.'//idchar,action='WRITE',&
  status='REPLACE', DELIM='QUOTE', iostat=ios)
  if (ios/=0) write(*, *) 'makerestartpoint: Warning: Failed opening', &
  parameterfile
  pwriter = new_parameter_writer(parameterunit)
  call mce_writeparameters(pwriter)
  close(parameterunit)

  !! Write configurations to a restartfile
  configurationunit = fileunit_getfreeunit()
  configurationfile='restartconfiguration.'//idchar
  open(file=configurationfile, unit=configurationunit, &
  action='WRITE', status='REPLACE', iostat=ios)
  if (ios/=0) write(*, *) 'makerestartpoint: Warning: Failed opening', &
  configurationfile

  call get_system(simbox, particles)
  call writestate(restartwriter, configurationunit, simbox, particles)
  close(configurationunit)
  deallocate(particles)

  !! Write random number generator state to a file.
  rngunit = fileunit_getfreeunit()
  rngfile='restartmtstate.'//trim(adjustl(idchar))
  open(unit=rngunit, file=rngfile, action='WRITE', status='REPLACE',&
  form='UNFORMATTED', iostat=ios)
  if(0 /= ios) then
    write(*, *) 'makerestartpoint: Warning: Failed opening ', rngfile
  else
    !call mtsave(rngunit, 'f')
    call save(mts(0), rngunit)
    close(rngunit)
  end if
end subroutine


!> All actions performed only during equilibration sweeps and not during
!! production sweeps should be gathered inside this routine.
!!  
subroutine runequilibrationtasks
  real(dp) :: temperature
  if (moveadjustperiod /= 0) then
    if (mod(isweep, moveadjustperiod) == 0) then
      call updatemaxvalues
    end if
  end if
  if (ptadjustperiod /= 0) then
    if (mod(isweep, ptadjustperiod) == 0) then
      temperature = gettemperature()
      call pt_adjusttemperature(temperature)
      call settemperature(temperature)
    end if
  end if
end subroutine 

!> All actions done during both the equilibration and production sweeps should
!! be gathered inside this routine for clarity
!! 
subroutine runproductiontasks
  integer :: ios
  character(len=80) :: parameterfile
  type(parameter_writer) :: writer
  type(factory) :: coordinatewriter
  type(poly_box) :: simbox
  type(particledat), allocatable :: particles(:)
  parameterfile='parameters.'//trim(adjustl(idchar))
  if (mod(isweep, productionperiod) == 0) then
    call get_system(simbox, particles)
    call writestate(coordinatewriter, coordinateunit, simbox, particles)
    deallocate(particles) 
    pwunit = fileunit_getfreeunit()
    open(UNIT=pwunit, FILE=parameterfile, action='WRITE', position='APPEND',&
    DELIM='QUOTE', IOSTAT=ios) 
    if (ios/=0) then
      write(*, *) 'runproductiontasks: error opening', parameterfile
      stop
    end if
    writer = new_parameter_writer(pwunit)
    call mce_writeparameters(writer)
    call delete(writer)
    close(pwunit)
    call resetcounters
  end if
end subroutine
  
end module
