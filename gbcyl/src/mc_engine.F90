!> A module for streering the simulation and controlling it's input and output.
module mc_engine
use nrtype
use utils
!use class_cylformatter
use class_factory
!use particlewall, only: initptwall, particlewall_writeparameters
!use gayberne, only: gayberne_init, gb_writeparameters
!use mtmod, only: sgrnd, mtsave, mtget, mts
use mt_stream
use particle, only: particledat, initParticle, particle_writeparameters
use class_poly_box
use openmpsweep, only: ompsweep => sweep, ompsweep_init
use mc_sweep, only: mc_sweep_init => init, updatemaxvalues, sweep, &
mc_sweep_writeparameters, resetcounters, gettemperature, settemperature 
!, getpressure
!use class_poly_nbrlist
!use energy
use pt
use m_fileunit
use class_parameterizer
use class_parameter_writer
implicit none
private

public :: init
public :: run
public :: finalize
public :: mce_writeparameters

interface init
  module procedure mce_init
end interface
  
type(particledat), dimension(:), pointer, save :: particles
type(poly_box), save :: simbox
integer, save :: nequilibrationsweeps = 0
integer, save :: nproductionsweeps = 0
integer, save :: productionperiod = 1
integer, save :: moveadjustperiod = 100
integer, save :: ptadjustperiod = 100
integer, save :: restartperiod = 10000
integer, save :: isweep = 0
integer, save :: rngunit
integer, save :: pwunit
character(len = 80), save :: rngfile = 'mtstate.'
character(len = 9), save :: idchar
type(factory), save :: coordinatereader
integer, save :: coordinateunit
logical :: isopenmp = .false.
type(mt_state), save :: mts
type(mt_state), save :: mts_new
character(len=80), save :: parameterfilename

contains
  
!> Initializes the state of the simulation. 
!!
!! pre-condition: state must not be initialized by other means. 
!! pre-condition: MPI_INIT has been called.
!! post-condition: simulation can be run.  
subroutine mce_init(id, ispt)
  integer, intent(in) :: id 
  logical, intent(in) :: ispt
  logical :: isrestart 
  character(len = 50) :: statefile 
  integer :: seed
  type(parameterizer) :: parameterreader
  integer :: ios
  type(parameter_writer) :: writer
  character(len=80) :: parameterinputfile
  write(idchar, '(I9)') id 
  idchar = trim(adjustl(idchar))
  isrestart = .false.
  parameterinputfile='restartparameters.'//trim(adjustl(idchar))
  parameterreader = new_parameterizer(parameterinputfile, ios=ios)
  if (ios/=0) then
    write(*, *) 'Could not open ', parameterinputfile
    parameterinputfile='parameters.in.'//trim(adjustl(idchar))
    parameterreader = new_parameterizer(parameterinputfile, ios=ios)
    if (ios/=0) then
      write(*,*) 'Could not open ', parameterinputfile, '. Stopping.'
      stop
    end if
  end if
  rngfile = trim(adjustl(rngfile)) // trim(adjustl(idchar))
  rngunit = fileunit_getfreeunit()
  !! For mtmod use
  !!open(unit = rngunit, file = trim(adjustl(rngfile)), &
  !!action = 'READ', status = 'OLD', form = 'FORMATTED', iostat = ios)
  !! For mt_stream use
  open(unit = rngunit, file = trim(adjustl(rngfile)), &
  action = 'READ', status = 'OLD', form = 'UNFORMATTED', iostat = ios)
  call set_mt19937
  call new(mts)
  if(0 /= ios) then
    write(*, *) 'Warning: Could not open ' // trim(rngfile) // &
    ' for reading!'
    call getparameter(parameterreader, 'seed', seed)
    !call sgrnd(seed + id)
    call init(mts, seed)
    if (id > 0) then 
      call create_stream(mts, mts_new, id)
      mts = mts_new
    end if
    !write(*, *) 'first random is:', genrand_double1(mts)
  else
    !call mtget(rngunit, 'f')
    call read(mts, rngunit)
    close(rngunit)
  end if
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
  !coordinatereader = new_cylformatter(statefile)
  call getparameter(parameterreader, 'isopenmp', isopenmp)
  call getparameter(parameterreader, 'restartperiod', restartperiod)
  !! Find last configuration written in the statefile.
  !call findlast(coordinatereader, isfound)
  !if (.not. isfound) then
  !  write(*, *) 'Error: Could not find a configuration from ' // &
  !  trim(statefile)
  !end if
  coordinateunit = fileunit_getfreeunit()
  statefile = 'restartconfiguration.'//trim(adjustl(idchar))
  open(file=statefile, unit=coordinateunit, action='READWRITE', iostat=ios)
  if (ios /= 0) then 
    write(*, *) statefile, ' not found'
    statefile = 'configurations.' // trim(adjustl(idchar))
    open(file=statefile, unit=coordinateunit, action='READWRITE', iostat=ios)
    if (ios /= 0) stop 'Could not open coordinate file'
  end if
  call readstate(coordinatereader, coordinateunit, simbox, particles, ios)
  if (ios > 0) then
    write(*, *) 'Error ', ios,' reading', statefile 
    stop
  end if
  !! Initialize modules. 
  if (isopenmp) then 
    call ompsweep_init(parameterreader, simbox, particles)
  else
    call mc_sweep_init(parameterreader, simbox, particles)
  end if
  call initparticle(parameterreader)
  if (ispt) call pt_init(parameterreader)
  call delete(parameterreader)

  !! Write initial parameters to parameters.out.(id)
  !! replacing old parameters.out.(id) file.
  pwunit = fileunit_getfreeunit()
  parameterfilename = 'parameters.out.' // idchar
  open(UNIT=pwunit, FILE = parameterfilename, action = 'WRITE', IOSTAT = ios)
  if (ios/=0) then 
    write(*, *) 'mce_init: error opening', parameterfilename
    stop
  end if
  writer = new_parameter_writer(pwunit)
  call mce_writeparameters(writer)
  call delete(writer)
  close(pwunit)
end subroutine 

!> Writes the parameters and observables of this module and its children
!! with the class_parameter_writer module and by calling the write routines
!! of the child modules.
!!
subroutine mce_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  integer :: ios
  call writecomment(writer, 'mc engine parameters')
  call writeparameter(writer, 'n_equilibration_sweeps', &
  nequilibrationsweeps)
  call writeparameter(writer, 'n_production_sweeps', nproductionsweeps)
  call writeparameter(writer, 'i_sweep', isweep)
  call writeparameter(writer, 'production_period', productionperiod)
  call writeparameter(writer, 'move_adjusting_period', moveadjustperiod)
  call writeparameter(writer, 'pt_adjusting_period', ptadjustperiod)
  call writeparameter(writer, 'isopenmp', isopenmp)
  call writeparameter(writer, 'restartperiod', restartperiod)
  rngunit = fileunit_getfreeunit()
  open(unit = rngunit, file = rngfile, action = 'WRITE', &
  status = 'REPLACE', form = 'UNFORMATTED', iostat = ios)
  if(0 /= ios) then
    write(*, *) 'Warning: Could not open rng file for writing!'
  else
    !call mtsave(rngunit, 'f')
    call save(mts, rngunit)
  end if
  !if (isopenmp) then
  !  call ompsweep_writeparameters(writer)
  !else
    call mc_sweep_writeparameters(writer)
  !end if
  call particle_writeparameters(writer)
  call pt_writeparameters(writer)
  close(rngunit)
end subroutine

!> Finalizes the simulation.
!!
!! pre-condition: simulation must be initialized.
!! post-conditions: 1. memory allocated by this program is freed.
!!                  2. program ends.
!!
subroutine finalize(id)
  integer, intent(in) :: id
  close(coordinateunit)
  call makerestartpoint
  if (associated(particles)) deallocate(particles)
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
    if (isopenmp) then 
      call ompsweep(simbox, particles, mts)
    else
      call sweep(simbox, particles, mts)
    end if
    if (isweep <= nequilibrationsweeps) then
      call runequilibrationtasks
    end if
    call runproductiontasks
    if (mod(isweep, restartperiod)==0) then
      call makerestartpoint()
    end if
  end do
end subroutine 

subroutine makerestartpoint
  type(factory) :: restartwriter
  integer :: configurationunit
  integer :: parameterunit
  type(parameter_writer) :: pwriter
  character(len=80) :: parameterfilename
  character(len=80) :: configurationfilename
  integer :: ios

  !! Write parameters to a restartfile
  parameterunit = fileunit_getfreeunit()
  open(UNIT=pwunit, FILE ='restartparameters.'//idchar, action='WRITE', &
  iostat=ios)
  if (ios/=0) write(*, *) 'makerestartpoint: Warning: Failed opening', &
  parameterfilename
  pwriter = new_parameter_writer(parameterunit)
  call mce_writeparameters(pwriter)
  close(parameterunit)

  !! Write configurations to a restartfile
  configurationunit = fileunit_getfreeunit()
  configurationfilename='restartconfiguration.'//idchar
  open(file=configurationfilename, unit=configurationunit, &
  action='WRITE', iostat=ios)
  if (ios/=0) write(*, *) 'makerestartpoint: Warning: Failed opening', &
  configurationfilename
  call writestate(restartwriter, configurationunit, simbox, particles)
  close(configurationunit)
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
  type(parameter_writer) :: writer
  if (mod(isweep, productionperiod) == 0) then
    call writestate(coordinatereader, coordinateunit, simbox, particles) 
    pwunit = fileunit_getfreeunit()
    open(UNIT = pwunit, FILE = parameterfilename, action = 'WRITE', & 
    position = 'APPEND', IOSTAT = ios) 
    if (ios/=0) then
      write(*, *) 'runproductiontasks: error opening', parameterfilename
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
