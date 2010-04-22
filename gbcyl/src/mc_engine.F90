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
getpressure, gettemperature, settemperature, mc_sweep_writeparameters
!!use verlet, only: initvlist, freevlist, verlet_writeparameters
use class_poly_nbrlist
use cell, only: list, writetostdout
use cell_energy, only: cell_energy_init, new_celllist, &
cell_energy_writeparameters => writeparameters
use energy, only: totalenergy
use pt
use mpi
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
integer, save :: adjustingperiod = 1
integer, save :: isweep = 0
!integer, save :: restartperiod = 1
integer, save :: energyunit
integer, save :: rngunit
integer, save :: pwunit
character(len = 80), save :: rngfile = 'mt-state.'
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
  rngunit = freeunit()
  write(idchar, '(I9)') id 
  idchar = trim(adjustl(idchar))
  isrestart = .false.
  do i = 0, ntasks - 1
    if ( i == id) then
      parameterreader = new_parameterizer('fin-kgb.' // &
      trim(adjustl(idchar)) // '.pms')
      rngfile = trim(adjustl(rngfile)) // trim(adjustl(idchar))
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
      !call getparameter(parameterreader, 'initstate', statefile)  
      statefile = 'fin-kgb.' // trim(adjustl(idchar)) // '.mol'
      call getparameter(parameterreader, 'n_equilibration_sweeps', &
      nequilibrationsweeps)
      call getparameter(parameterreader, 'n_production_sweeps', &
      nproductionsweeps)
      call getparameter(parameterreader, 'production_period', &
      productionperiod)
      call getparameter(parameterreader, 'i_sweep', isweep)
      call getparameter(parameterreader, 'adjusting_period', adjustingperiod)
      cf = new_cylformatter(statefile)
      !call findlast(cf, isfound)
      !if (.not. isfound) then
      !  write(*, *) 'Error: Could not find a configuration from ' // &
      !  trim(statefile)
      !end if
      call readstate(cf, particles, nparticles, radius, height)
      !! Initialize modules. 
      simbox = new_cylinder(2._dp * radius, height)
      call initptwall(parameterreader)
      call gayberne_init(parameterreader)
      call cell_energy_init(parameterreader)
      nbrlist = new_celllist(simbox, particles(1:nparticles))
      call mc_sweep_init(parameterreader)
      !call initvlist(particles, nparticles, simbox, parameterreader)
      call initparticle(parameterreader)
      call openenergyfile
      call pt_init()
      call totalenergy(simbox, particles(1:nparticles), nbrlist, etot, overlap)
      if (overlap) then
        stop 'Overlap found in starting configuration. Simulation will stop.'
      else
        write(*, *) 'Total energy of initial configuration is ', etot
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, rc)
      call delete(parameterreader)
    end if
  end do
end subroutine 

subroutine openenergyfile
  integer :: filestatus
  character(len = 30) :: energyfile
  energyfile = 'thermodynamics.dat.' // idchar
  energyunit = freeunit()
  open(unit = energyunit, file = energyfile, action = 'WRITE', &
   & status = 'REPLACE', delim='QUOTE', iostat = filestatus)
  if (filestatus /= 0) then
    stop 'Could not open energy file for writing! Stopping.'
  end if
end subroutine 

subroutine writerestart
  type(parameter_writer) :: writer
  character(len = 80) :: parameterfilename
  integer :: ios
  pwunit = freeunit()
  parameterfilename = 'parameters.out.' // idchar
  open(UNIT = pwunit, FILE = parameterfilename, action = 'WRITE', & 
  status = 'REPLACE', IOSTAT = ios) 
  writer = new_parameter_writer(pwunit)
  call writecomment(writer, 'mc engine parameters')
  call writeparameter(writer, 'is_restart', .true.)
  call writeparameter(writer, 'initstate', '\"simdata.out.' // idchar &
  //'\"')  
  call writeparameter(writer, 'n_equilibration_sweeps', &
  nequilibrationsweeps)
  call writeparameter(writer, 'n_production_sweeps', nproductionsweeps)
  call writeparameter(writer, 'i_sweep', isweep)
  call writeparameter(writer, 'production_period', productionperiod)
  open(unit = rngunit, file = rngfile, action = 'WRITE', &
  status = 'REPLACE', form = 'FORMATTED', iostat = ios)
  if(0 /= ios) then
    write(*, *) 'Warning: Could not open rng file for writing!'
  else
    call mtsave(rngunit, 'f')
  end if
  !! Initialize modules. 
  call particlewall_writeparameters(writer)
  call gb_writeparameters(writer)
  call mc_sweep_writeparameters(writer)
!  call verlet_writeparameters(writer)
  call cell_energy_writeparameters(writer)
  call particle_writeparameters(writer)
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
!  call freevlist()
  if (associated(particles)) deallocate(particles)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, rc)
  if (myid == 0) write (*, *) 'End of program gbcyl.'
  close(energyunit)
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
    if (isweep .le. nequilibrationsweeps) then
      call runequilibrationtasks
    end if
    call runproductiontasks
  end do
end subroutine 
  
subroutine runequilibrationtasks
  real(dp) :: temperature
  if (mod(isweep, adjustingperiod) == 0) then
    !call updatemaxvalues(nparticles, adjustingperiod)
    temperature = gettemperature()
    call pt_adjusttemperature(temperature)
    call settemperature(temperature)
  end if
end subroutine 

subroutine runproductiontasks
  if (mod(isweep, productionperiod) == 0) then
    !! :TODO: Make i/o compatible to poly_box   
    call writestate(cf, particles, nparticles, getx(simbox)/2._dp, &
    getz(simbox))
    call writethermodynamics
    call writerestart
  end if
end subroutine
  
subroutine writethermodynamics
  real(dp) :: etot
  logical :: overlap
  call totalenergy(simbox, particles(1:nparticles), nbrlist, etot, overlap)
  if (overlap) then 
    stop 'Writing statistics with an overlapping configuration!'
  end if
  write(energyunit, '(' // fmt_char_int() // ', 3' // fmt_char_dp() // &
  ')') isweep, etot, volume(simbox), etot + getpressure() * &
  volume(simbox) 
end subroutine 

end module mc_engine
