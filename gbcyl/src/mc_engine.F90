module mc_engine
  use nrtype, only: dp
  use io, only: ReadParams, readstate, writestate, init_io => init, &
    finalizeOutput, io_save_state => save_state, io_load_state => load_state
  use particlewall, only: initptwall, particlewall_save_state => save_state, &
    particlewall_load_state => load_state
  use gayberne, only: gayberne_init => init, &
    gayberne_save_state => save_state, gayberne_load_state => load_state
  use mtmod, only: sgrnd
  use particle, only: particledat, initParticle, &
    particle_save_state => save_state, &
    particle_load_state => load_state
  use cylinder, only: initcylinder, getHeight, getRadius, &
    cylinder_save_state => save_state, &
    cylinder_load_state => load_state
  use mc_sweep, only: mc_sweep_init => init, updatemaxvalues, sweep, &
    mc_sweep_save_state => save_state, &
    mc_sweep_load_state => load_state
  use verlet, only: initvlist, freevlist, &
    verlet_save_state => save_state, &
    verlet_load_state => load_state
  implicit none


  public :: init
  public :: run
  public :: finalize
  public :: write_restart
  public :: read_restart
  public :: write_restart_to 
  public :: read_restart_from
  public :: equilibration_sweeps
  public :: production_sweeps  



  private
  
  integer, save :: n_particles_
  type(particledat), dimension(:), pointer, save :: particles_
  integer, save :: n_equilibration_sweeps_
  integer, save :: n_production_sweeps_
  integer, save :: production_period_
  integer, save :: adjusting_period_
  integer, save :: i_sweep_
  integer, save :: restart_period_
  integer, save :: restart_unit_
  character(len=*), parameter :: restart_file_ = "restart.gbcyl"
  namelist /mc_engine_nml/ n_particles_, n_equilibration_sweeps_, &
    n_production_sweeps_, production_period_, i_sweep_, restart_period_, &
    restart_unit_, adjusting_period_ 

  

  contains

  !! Initializes the state of the simulation. 
  !!
  !! pre-condition: state must not be initialized by other means. 
  !! post-condition: simulation can be run. 
  !! 
  subroutine init
    implicit none
    character(len = 78) :: statefile 
    real(dp) :: temperature
    real(dp) :: pressure
    integer :: anchor
    integer :: volume_scaling_type
    real(dp) :: Kw
    integer :: seed
    real(dp) :: kappa_epsilon_sgb
    real(dp) :: epsilon_0_sgb
    real(dp) :: r_sphere
    real(dp) :: spmyy
    real(dp) :: epsilon_ss
    real(dp) :: sigma_0_sgb
    real(dp) :: kappa_sigma_sgb
    integer :: allign
    integer :: debug 
    real(dp) :: radius
    real(dp) :: height
    real(dp), parameter :: maxtrans = 0.156
    real(dp), parameter :: maxangle = 0.170
    real(dp), parameter :: kappa_sigma = 4.4
    real(dp), parameter :: kappa_epsilon = 20.0
    real(dp), parameter :: mu = 1.0
    real(dp), parameter :: nu = 1.0
    real(dp), parameter :: sigma_0 = 1.0
    real(dp), parameter :: epsilon_0 = 1.0
    real(dp), dimension(3), parameter :: magnetic_field_direction = &
      & (/0.0, 0.0, 1.0/)
    real(dp), parameter :: magnetic_field_teslas = 11.74
    logical, parameter :: is_magnet_on = .false.
    call ReadParams(statefile, n_equilibration_sweeps_, n_production_sweeps_, &
      & production_period_, temperature, pressure, anchor, & 
      & volume_scaling_type, Kw, seed, kappa_epsilon_sgb, epsilon_0_sgb, & 
      & r_sphere, spmyy, epsilon_ss, sigma_0_sgb, kappa_sigma_sgb, allign, & 
      & debug)
    call readstate(statefile, particles_, n_particles_, radius, height); 
    !! Initialize modules. 
    call initptwall(anchor, Kw)
    call sgrnd(seed)  
    call gayberne_init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
    call init_io
    call initcylinder(radius, height)
    call mc_sweep_init(volume_scaling_type, temperature, pressure)
    call initvlist(particles_, n_particles_)
    call initParticle(maxtrans, maxangle)
    i_sweep_ = 0
    restart_period_ = 1000
    restart_unit_ = 13
    adjusting_period_ = 20
  end subroutine init



  !! Finalizes the simulation.
  !!
  !! pre-condition: simulation must be initialized.
  !! post-conditions: 1. memory allocated by this program is freed.
  !!                  2. program ends.
  !!
  subroutine finalize()
    implicit none
    call finalizeOutput
    call freevlist()
    if (associated(particles_)) deallocate(particles_)
    write (*, *) 'End of program gbcyl.'
  end subroutine finalize


  
  subroutine equilibration_sweeps(n_equilibration_sweeps)
    implicit none
    integer, intent(in) :: n_equilibration_sweeps
    n_equilibration_sweeps_ = n_equilibration_sweeps
  end subroutine equilibration_sweeps



  subroutine production_sweeps(n_production_sweeps)
    implicit none
    integer, intent(in) :: n_production_sweeps
    n_production_sweeps_ = n_production_sweeps
  end subroutine production_sweeps



  !! Writes the state of the simulation to the restart_file_. 
  !! 
  !! :NOTE: This operation has to be symmetric with the operation 
  !! read_restart, so that the modules are saved in the same order as they 
  !! are loaded.  
  !!
  !! pre-conditions: 1. simulation must be initialized.
  !!                 2. restart_file_ must not be in use. 
  !! 
  subroutine write_restart
    implicit none
    integer :: file_status
    open(unit = restart_unit_, file = restart_file_, action = 'WRITE', &
     & status = 'REPLACE', iostat = file_status)
    if (file_status /= 0) then
      stop 'Could not open restart file for writing! Stopping.'
    end if
    call write_restart_to(restart_unit_)
    close(restart_unit_)
  end subroutine write_restart



  !! Writes the state of the simulation to the @p write_unit. 
  !! 
  !! :NOTE: This operation has to be symmetric with the operation load_state, 
  !! so that the modules are saved in the same order as they are loaded.
  !!
  !! pre-conditions: 1. @p write_unit must be ready for writing.
  !!                 2. state of simulation must be initialized.
  !!
  !! @param write_unit the writing unit where the state will be saved.
  !!
  subroutine write_restart_to(write_unit)
    implicit none
    integer, intent(in) :: write_unit
    call io_save_state(write_unit)
    call mc_sweep_save_state(write_unit)
    call cylinder_save_state(write_unit)
    call particlewall_save_state(write_unit)
    call gayberne_save_state(write_unit)
    call particle_save_state(write_unit)
    call save_state(write_unit)
    call verlet_save_state(write_unit)
  end subroutine write_restart_to



  !! Reads The state of the simulation from the restart_file_.
  !! 
  !! pre-condition: 1. restart_file_ must not be in use.
  !!                2. simulation state must not be initialized.
  !! 
  !! @see read_restart_from(read_unit)
  !!
  subroutine read_restart
    implicit none
    integer :: file_status
    open(unit = restart_unit_, file = restart_file_, action = 'READ', &
      & status = 'OLD', iostat = file_status)
    if (file_status /= 0) then
      stop 'Could not open restart file for reading! Stopping.'
    end if
    call read_restart_from(restart_unit_)
    close(restart_unit_)
  end subroutine read_restart

  

  !! Reads The state of the simulation from @p read_unit.
  !!
  !! :NOTE: This operation has to be symmetric with the operation 
  !! write_restart_to, so that the state of the modules is loaded in the same 
  !! order as saved.
  !!
  !! pre-condition: @p read_unit must be ready for reading.
  !!
  !! @p read_unit the reading unit where the state is read from.
  !!
  subroutine read_restart_from(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    call io_load_state(read_unit)
    call mc_sweep_load_state(read_unit)
    call cylinder_load_state(read_unit)
    call particlewall_load_state(read_unit)
    call gayberne_load_state(read_unit)
    call particle_load_state(read_unit)
    call load_state(read_unit)
    call verlet_load_state(read_unit)
  end subroutine read_restart_from



  !! Saves the state of this module to @p write_unit
  !! 
  !! @param write_unit the unit to write the state to.
  !!
  subroutine save_state(write_unit)
    implicit none
    integer, intent(in) :: write_unit
    write(write_unit, NML = mc_engine_nml) 
    write(write_unit, *) particles_(1:n_particles_) 
  end subroutine save_state



  !! Loads the state of this module from @p read_unit.
  !!
  !! @param read_unit the unit to read the state from.
  !!
  subroutine load_state(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    read(read_unit, NML = mc_engine_nml)
    if(.not. associated(particles_)) allocate(particles_(n_particles_))
    read(read_unit, *) particles_(1:n_particles_)
  end subroutine load_state



  !! Runs the simulation. 
  !! 
  !! pre-condition: simulation state has to be initialized. 
  !! 
  subroutine run
    implicit none
    do while ((n_equilibration_sweeps_ + n_production_sweeps_) > i_sweep_)
      if (mod(i_sweep_, restart_period_) == 0) call write_restart
      if (i_sweep_ .lt. n_equilibration_sweeps_) call run_equilibration_tasks
      if (i_sweep_ .ge. n_equilibration_sweeps_) call run_production_tasks
      i_sweep_ = i_sweep_ + 1
      call sweep(particles_, n_particles_)
    end do
  end subroutine run


  
  subroutine run_equilibration_tasks
    implicit none 
    if (mod(i_sweep_, adjusting_period_) == 0) then
      call updateMaxValues(n_particles_, adjusting_period_)
    end if
  end subroutine run_equilibration_tasks



  subroutine run_production_tasks
    implicit none
    real(dp) :: radius
    real(dp) :: height
    if (mod(i_sweep_, production_period_) == 0) then
      radius = getRadius()
      height = getHeight()
      call writestate(particles_, n_particles_, radius, height)
    end if
  end subroutine run_production_tasks



end module mc_engine
