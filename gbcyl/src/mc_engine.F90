module mc_engine
  use nrtype
  use utils
  use io, only: readstate, writestate, init_io => init, finalizeOutput, &
    io_write_parameters
  use particlewall, only: initptwall, particlewall_write_parameters
  use gayberne, only: gayberne_init => init, gb_write_parameters
  use mtmod, only: sgrnd, mtsave, mtget
  use particle, only: particledat, initParticle, particle_write_parameters
  use class_poly_box
  use cylinder, only: new_cylinder
  use mc_sweep, only: mc_sweep_init => init, updatemaxvalues, sweep, pressure, &
    mc_sweep_write_parameters
  use verlet, only: initvlist, freevlist, verlet_write_parameters
  use energy, only: total_energy
  use pt
  use mpi
  use class_parameterizer
  use class_parameter_writer
  implicit none

  public :: init
  public :: run
  public :: finalize
  public :: write_restart

  private
  
  integer, save :: n_particles_
  type(particledat), dimension(:), pointer, save :: particles_
  type(poly_box), save :: simbox_
  integer, save :: n_equilibration_sweeps_ = 0
  integer, save :: n_production_sweeps_ = 0
  integer, save :: production_period_ = 1
  integer, save :: adjusting_period_ = 1
  integer, save :: i_sweep_ = 0
  integer, save :: restart_period_ = 1
  integer, parameter :: energy_unit_ = 15
  integer, parameter :: rng_unit_ = 10
  character(len = *), parameter :: rng_file_ = 'mt-restart.dat'
  character(len = 9), save :: id_char_

  contains

  !! Initializes the state of the simulation. 
  !!
  !! pre-condition: state must not be initialized by other means. 
  !! post-condition: simulation can be run. 
  !! 
  subroutine init
    implicit none
    logical :: is_restart 
    character(len = 50) :: statefile 
    real(dp) :: radius
    real(dp) :: height
    logical :: overlap
    real(dp) :: e_tot
    integer :: id 
    integer :: rc
    integer :: n_tasks
    integer :: i
    integer :: seed
    type(parameterizer) :: parameter_reader
    integer :: ios
    call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, rc)
    write(id_char_, '(I9)') id 
    id_char_ = trim(adjustl(id_char_))
    is_restart = .false.
    do i = 0, n_tasks - 1
      if ( i == id) then
        parameter_reader = new_parameterizer('fin-kgb.' // trim(adjustl(id_char_)) // '.pms')
        open(unit = rng_unit_, file = rng_file_, action = 'READ', &
        status = 'OLD', form = 'FORMATTED', iostat = ios)
        if(0 /= ios) then
          write(*, *) 'Warning: Could not open rng file for reading!'
          call get_parameter(parameter_reader, 'seed', seed)
          call sgrnd(seed)
        else
          call mtget(rng_unit_, 'f')
          close(rng_unit_)
        end if
        !call get_parameter(parameter_reader, 'initstate', statefile)  
        statefile = 'fin-kgb.' // trim(adjustl(id_char_)) // '.mol'
        call get_parameter(parameter_reader, 'n_equilibration_sweeps', n_equilibration_sweeps_)
        call get_parameter(parameter_reader, 'n_production_sweeps', n_production_sweeps_)
        call get_parameter(parameter_reader, 'production_period', production_period_)
        call get_parameter(parameter_reader, 'i_sweep', i_sweep_)
        call readstate(statefile, particles_, n_particles_, radius, height)
        !! Initialize modules. 
        call initptwall(parameter_reader)
        call gayberne_init(parameter_reader)
        call init_io(statefile)
        simbox_ = new_cylinder(2._dp * radius, height)
        call mc_sweep_init(parameter_reader)
        call initvlist(particles_, n_particles_, simbox_, parameter_reader)
        call initparticle(parameter_reader)
        restart_period_ = 1000
        adjusting_period_ = 20
        call open_energyfile
        call pt_init()
        call total_energy(particles_, n_particles_, simbox_, e_tot, overlap)
        if (overlap) then
          stop 'Overlap found in starting configuration. Simulation will stop.'
        else
          write(*,*) 'Total energy of initial configuration is ', e_tot
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, rc)
        call delete(parameter_reader)
      end if
    end do
  end subroutine 

  subroutine open_energyfile
    implicit none
    integer :: file_status
    character(len = 30) :: energy_file
    energy_file = 'thermodynamics.dat.' // trim(adjustl(id_char_))
    open(unit = energy_unit_, file = energy_file, action = 'WRITE', &
     & status = 'REPLACE', delim='QUOTE', iostat = file_status)
    if (file_status /= 0) then
      stop 'Could not open energy file for writing! Stopping.'
    end if
  end subroutine 

  subroutine write_restart
    type(parameter_writer) :: writer
    integer, parameter :: pw_unit = 11
    character(len = 80) :: parameter_filename
    integer :: ios
    parameter_filename = 'parameters.out.' // id_char_
    open(UNIT = pw_unit, FILE = parameter_filename, action = 'WRITE', & 
    status = 'REPLACE', IOSTAT = ios) 
    writer = new_parameter_writer(pw_unit)
    call write_comment(writer, 'mc engine parameters')
    call write_parameter(writer, 'is_restart', .true.)
    call write_parameter(writer, 'initstate', '\"simdata.out.' // id_char_ //'\"')  
    call write_parameter(writer, 'n_equilibration_sweeps', n_equilibration_sweeps_)
    call write_parameter(writer, 'n_production_sweeps', n_production_sweeps_)
    call write_parameter(writer, 'i_sweep', i_sweep_)
    call write_parameter(writer, 'production_period', production_period_)
!    call write_state(state_writer, particles_, n_particles_, simbox_)
    open(unit = rng_unit_, file = rng_file_, action = 'WRITE', &
    status = 'REPLACE', form = 'FORMATTED', iostat = ios)
    if(0 /= ios) then
      write(*, *) 'Warning: Could not open rng file for writing!'
    else
      call mtsave(rng_unit_, 'f')
    end if
    !! Initialize modules. 
    call particlewall_write_parameters(writer)
    call gb_write_parameters(writer)
    !call write_io How to write io state?
    !simbox_ = new_cylinder(2._dp * radius, height) How to handle this when 
    ! writing? How about making a new geometry type independent coordinate file?
    call mc_sweep_write_parameters(writer)
    call verlet_write_parameters(writer)
    call particle_write_parameters(writer)
    call delete(writer)
    close(rng_unit_)
  end subroutine

  !! Finalizes the simulation.
  !!
  !! pre-condition: simulation must be initialized.
  !! post-conditions: 1. memory allocated by this program is freed.
  !!                  2. program ends.
  !!
  subroutine finalize()
    integer :: my_id, rc
    call finalizeOutput
    call freevlist()
    if (associated(particles_)) deallocate(particles_)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, rc)
    if (my_id == 0) write (*, *) 'End of program gbcyl.'
    close(energy_unit_)
    call pt_finalize()
    !close(molecule_unit_)
  end subroutine 
  
  !! Runs the simulation. 
  !! 
  !! pre-condition: simulation state has to be initialized. 
  !! 
  subroutine run
    call write_restart
    do while (i_sweep_ < n_equilibration_sweeps_ + n_production_sweeps_)
      i_sweep_ = i_sweep_ + 1 
      call sweep(particles_, n_particles_, simbox_)
      if (i_sweep_ .le. n_equilibration_sweeps_) then
        call run_equilibration_tasks
      else 
        call run_production_tasks
      end if
      if (mod(i_sweep_, restart_period_) == 0) call write_restart
    end do
    call write_restart
  end subroutine 
  
  subroutine run_equilibration_tasks
    if (mod(i_sweep_, adjusting_period_) == 0) then
      call updateMaxValues(n_particles_, adjusting_period_)
    end if
    if (mod(i_sweep_, production_period_) == 0) call write_thermodynamics
  end subroutine 

  subroutine run_production_tasks
    if (mod(i_sweep_, production_period_) == 0) then
      call writestate(particles_, n_particles_, get_x(simbox_)/2._dp, & ! :TODO: Make i/o compatible to poly_box 
        get_z(simbox_))
      call write_thermodynamics
    end if
  end subroutine 
  
  subroutine write_thermodynamics
    real(dp) :: e_tot
    logical :: overlap
    call total_energy(particles_, n_particles_, simbox_, e_tot, overlap)
    if (overlap) then
      stop 'Writing statistics with an overlapping configuration!'
    end if
    write(energy_unit_, '(' // fmt_char_int() // ', 3' // fmt_char_dp() // ')') i_sweep_, e_tot, volume(simbox_), &
      e_tot + pressure() * volume(simbox_) 
    write(*, '(' // fmt_char_int() // ', 3' // fmt_char_dp() // ')') i_sweep_, e_tot, volume(simbox_), &
      e_tot + pressure() * volume(simbox_) 
    !write(*, '(' // fmt_char_dp() // ')') i_sweep_, e_tot, volume(simbox_), &
    !  e_tot + pressure() * volume(simbox_) 
  end subroutine 


end module mc_engine
