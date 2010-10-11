module mc_engine
use nrtype
use utils
use class_cylformatter
use particlewall, only: initptwall, particlewall_writeparameters
use gayberne, only: gayberne_init, gb_writeparameters
use mtmod, only: sgrnd, mtsave, mtget
use particle, only: particledat, initParticle, particle_writeparameters
use class_poly_box
use cylinder, only: new_cylinder
use mc_sweep, only: mc_sweep_init => init, updatemaxvalues, sweep, &
getpressure, gettemperature, settemperature, mc_sweep_writeparameters, &
resetcounters
use class_poly_nbrlist
use energy
use pt
use mpi
use m_fileunit
use class_parameterizer
use class_parameter_writer
implicit none
private

public :: init
public :: run
public :: finalize
public :: writerestart

interface init
  module procedure mce_init
end interface
  
integer, save :: nparticles
type(particledat), dimension(:), pointer, save :: particles
type(poly_box), save :: simbox
integer, save :: nequilibrationsweeps = 0
integer, save :: nproductionsweeps = 0
integer, save :: productionperiod = 1
integer, save :: moveadjustperiod = 100
integer, save :: ptadjustperiod = 100
integer, save :: isweep = 0
integer, save :: rngunit
integer, save :: pwunit
character(len = 80), save :: rngfile = 'mtstate.'
character(len = 9), save :: idchar
type(cylformatter), save :: cf
type(poly_nbrlist), save :: nbrlist

contains
  
!> Initializes the state of the simulation. 
!!
!! pre-condition: state must not be initialized by other means. 
!! post-condition: simulation can be run. 
!< 
subroutine mce_init
  logical :: isrestart 
  character(len = 50) :: statefile 
  real(dp) :: radius
  real(dp) :: height
  logical :: overlap
  real(dp) :: etot
  integer :: id 
  integer :: rc
  integer :: ntasks
  integer :: i
  integer :: seed
  type(parameterizer) :: parameterreader
  integer :: ios
  logical :: isfound 
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
  write(idchar, '(I9)') id 
  idchar = trim(adjustl(idchar))
  isrestart = .false.
  do i = 0, ntasks - 1
    if ( i == id) then
      parameterreader = new_parameterizer('parameters.in.' // &
      trim(adjustl(idchar)))
      rngfile = trim(adjustl(rngfile)) // trim(adjustl(idchar))
      rngunit = fileunit_getfreeunit()
      open(unit = rngunit, file = trim(adjustl(rngfile)), &
      action = 'READ', status = 'OLD', form = 'FORMATTED', iostat = ios)
      if(0 /= ios) then
        write(*, *) 'Warning: Could not open ' // trim(rngfile) // &
        ' for reading!'
        call getparameter(parameterreader, 'seed', seed)
        call sgrnd(seed + id)
      else
        call mtget(rngunit, 'f')
        close(rngunit)
      end if
      statefile = 'configurations.' // trim(adjustl(idchar))
      call getparameter(parameterreader, 'n_equilibration_sweeps', &
      nequilibrationsweeps)
      call getparameter(parameterreader, 'n_production_sweeps', &
      nproductionsweeps)
      call getparameter(parameterreader, 'production_period', &
      productionperiod)
      call getparameter(parameterreader, 'i_sweep', isweep)
      call getparameter(parameterreader, 'move_adjusting_period', moveadjustperiod)
      call getparameter(parameterreader, 'pt_adjusting_period', ptadjustperiod)
      cf = new_cylformatter(statefile)
      !call findlast(cf, isfound)
      !if (.not. isfound) then
      !  write(*, *) 'Error: Could not find a configuration from ' // &
      !  trim(statefile)
      !end if
      call readstate(cf, particles, nparticles, radius, height)
      !! Initialize modules. 
      simbox = new_cylinder(2._dp * radius, height)
      call setid(simbox, id)
      call initptwall(parameterreader)
      call gayberne_init(parameterreader)
      call poly_nbrlist_init(parameterreader)
      nbrlist = create_nbrlist(simbox, particles)
      call mc_sweep_init(parameterreader)
      call initparticle(parameterreader)
      call pt_init()
      call potentialenergy(simbox, particles(1:nparticles), nbrlist, etot, overlap)
      if (overlap) then
        stop 'Overlap found in starting configuration. Simulation will stop.'
      else
        write(*, *) 'Total energy of initial configuration is ', etot
      end if
      call delete(parameterreader)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, rc)
  end do
end subroutine 

subroutine writerestart
  type(parameter_writer) :: writer
  character(len = 80) :: parameterfilename
  integer :: ios
  pwunit = fileunit_getfreeunit()
  parameterfilename = 'parameters.out.' // idchar
  open(UNIT = pwunit, FILE = parameterfilename, action = 'WRITE', & 
  position = 'APPEND', IOSTAT = ios) 
  writer = new_parameter_writer(pwunit)
  call writecomment(writer, 'mc engine parameters')
  call writeparameter(writer, 'n_equilibration_sweeps', &
  nequilibrationsweeps)
  call writeparameter(writer, 'n_production_sweeps', nproductionsweeps)
  call writeparameter(writer, 'i_sweep', isweep)
  call writeparameter(writer, 'production_period', productionperiod)
  call writeparameter(writer, 'move_adjusting_period', moveadjustperiod)
  call writeparameter(writer, 'pt_adjusting_period', ptadjustperiod)
  rngunit = fileunit_getfreeunit()
  open(unit = rngunit, file = rngfile, action = 'WRITE', &
  status = 'REPLACE', form = 'FORMATTED', iostat = ios)
  if(0 /= ios) then
    write(*, *) 'Warning: Could not open rng file for writing!'
  else
    call mtsave(rngunit, 'f')
  end if
  call particlewall_writeparameters(writer)
  call gb_writeparameters(writer)
  call mc_sweep_writeparameters(writer)
  call poly_nbrlist_writeparameters(writer)
  call particle_writeparameters(writer)
  call pt_writeparameters(writer)
  call delete(writer)
  close(rngunit)
end subroutine

!! Finalizes the simulation.
!!
!! pre-condition: simulation must be initialized.
!! post-conditions: 1. memory allocated by this program is freed.
!!                  2. program ends.
!!
subroutine finalize()
  integer :: myid, rc
  call delete(cf)
  call delete(nbrlist)
  if (associated(particles)) deallocate(particles)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, rc)
  call MPI_BARRIER(MPI_COMM_WORLD, rc)
  if (myid == 0) write (*, *) 'End of program gbcyl.'
  call pt_finalize()
end subroutine 
  
!! Runs the simulation. 
!! 
!! pre-condition: simulation state has to be initialized. 
!! 
subroutine run
  do while (isweep < nequilibrationsweeps + nproductionsweeps)
    isweep = isweep + 1 
    call sweep(simbox, particles(1:nparticles), nbrlist)
    if (isweep <= nequilibrationsweeps) then
      call runequilibrationtasks
    end if
    call runproductiontasks
  end do
end subroutine 
  
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

subroutine runproductiontasks
  if (mod(isweep, productionperiod) == 0) then
    !! :TODO: Make i/o compatible to poly_box   
    call writestate(cf, particles, nparticles, getx(simbox)/2._dp, &
    getz(simbox), getid(simbox))
    call writerestart
    call resetcounters
  end if
end subroutine
  
end module
