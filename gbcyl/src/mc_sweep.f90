!> Implements Metropolis Monte Carlo updates of the coordinates of the
!! particles and the simulation box dimensions. The module also
!! controls the parallel tempering updates. If compiled and run
!! with OpenMP, the domain decomposition algorithm is used for the
!! particle moves and energy calculations. 
module mc_sweep
  use nrtype
  use energy
  use class_poly_box
  use particle
  use particle_mover
  use beta_exchange, beta_exchange_init => init, be_finalize => finalize
  use class_parameterizer
  use class_parameter_writer
  use genvoltrial
  use utils 
  !$ use omp_lib
  use class_simplelist
  include 'rng.inc'
  implicit none  
  private 

  public :: mcsweep_init
  public :: sweep
  public :: updatemaxvalues
  public :: getpressure
  public :: gettemperature
  public :: mc_sweep_writeparameters
  public :: settemperature
  public :: resetcounters
  public :: movevol
  public :: pt_period
  public :: set_system, get_system
  public :: get_total_energy
  public :: test_configuration
  public :: mcsweep_finalize

  !> The number of accepted particle moves
  integer, save :: nacceptedmoves

  !> The number of accepted volume moves
  integer, save :: nacceptedscalings

  !> The maximum absolute change of volume in a trial volume update.
  real(dp), save :: maxscaling = 100._dp

  !> The simulation temperature.
  real(dp), save :: temperature = -1._dp


  !> The simulation pressure. Only meaningful in a constant-pressure
  !! simulation.
  real(dp), save :: pressure = -1._dp

  !> The desired acceptance ratio for trial particle moves.
  real(dp), save :: moveratio

  !> The desired acceptance ratio for trial volume updates.
  real(dp), save :: scalingratio

  !> The number of sweeps between parallel tempering updates.
  integer, save :: pt_period = 1

  !> The total energy.
  real(dp), save :: etotal = 0._dp

  !> Current volume of the system.
  real(dp), save :: currentvolume = 0._dp

  !> Counter for trial particle moves.
  integer, save :: nmovetrials = 0

  !> Counter for trial volume updates.
  integer, save :: nscalingtrials = 0

  !> The cell list type of neighbourlist.
  type(simplelist), save :: sl

  !> The types of trial volume updates.
  character(len = 200), dimension(:), allocatable, save :: scalingtypes

  !> True if the module is correctly initialized.
  logical :: is_initialized = .false.

  !> The simulation box in which the particles reside.
  type(poly_box), save :: simbox

  !> The particles. 
  type(particledat), allocatable, save :: particles(:)

  contains

!> Initializes the module by getting the module parameters from the
!! @p reader object.
!!
!! @param reader the object which gets parameters by name e.g. from an
!!        input file.
!! @param the_simbox the simulation box defining the boundary
!!        conditions and volume.
!! @param the_particles the particles to be updated with the MC moves.
!!
subroutine mcsweep_init(reader, the_simbox, the_particles)
  type(parameterizer), intent(in) :: reader
  type(poly_box), intent(in) :: the_simbox
  type(particledat), dimension(:), intent(in) :: the_particles
  real(dp) :: min_cell_length
  character(len = 200), save :: scalingtype = "z"
  call particlemover_init(reader)
  call energy_init(reader)
  call getparameter(reader, 'scaling_type', scalingtype)
  call parsescalingtype(scalingtype)
  call getparameter(reader, 'temperature', temperature)
  if (temperature < 0._dp) then
     write(*, *) 'mc_sweep: mcsweep_init: trying to set a negative temperature, stopping.'
     stop  
  end if
  call beta_exchange_init(1._dp / temperature)
  call getparameter(reader, 'pressure', pressure)
  if (pressure < 0._dp) then
     write(*, *) 'mc_sweep: mcsweep_init: trying to set a negative pressure, stopping.'
     stop  
  end if
  call getparameter(reader, 'move_ratio', moveratio)
  call getparameter(reader, 'scaling_ratio', scalingratio)
  call getparameter(reader, 'max_scaling', maxscaling)
  call getparameter(reader, 'pt_period', pt_period)
  call getparameter(reader, 'nmovetrials', nmovetrials)
  call getparameter(reader, 'nscalingtrials', nscalingtrials)
  call getparameter(reader, 'nacceptedscalings', nacceptedscalings)
  !! The initializations below may affect restart so that it does not result
  !! in the same simulation. 
  nacceptedmoves = 0
  nacceptedscalings = 0
  min_cell_length = get_cutoff() + 2._dp * get_max_translation()
  !call simplelist_init(reader)
  !$ if (.true.) then
  !$ call new_simplelist(sl, the_simbox, the_particles, min_cell_length, &
  !$& is_x_even = isxperiodic(simbox), is_y_even = isyperiodic(simbox), &
  !$& is_z_even = iszperiodic(simbox), cutoff=get_cutoff())
  !$ else 
  call new_simplelist(sl, the_simbox, the_particles, min_cell_length, &
       cutoff=get_cutoff())
  !$ end if
  call set_system(the_simbox, the_particles)
  is_initialized = .true.
end subroutine mcsweep_init

!> Finalizes the module state.
subroutine mcsweep_finalize
  call be_finalize
  if (allocated(particles)) deallocate(particles)
  call simplelist_deallocate(sl)
  if (allocated(scalingtypes)) deallocate(scalingtypes)
  is_initialized = .false.
end subroutine mcsweep_finalize

!> Returns the total energy of the system.
real(dp) function get_total_energy()
  if (.not. is_initialized) then
     stop 'mc_sweep:get_total_energy: Error: module not initialized!'
  end if
  get_total_energy = etotal
end function get_total_energy

!> Checks and parses the @p scalingtype string to the scalingtypes
!! array. 
subroutine parsescalingtype(scalingtype)
  character(len=*), intent(in) :: scalingtype
  logical :: isok
  isok = (0 == verify(trim(adjustl(scalingtype)), 'xyz,') .and. &
       0 /= len_trim(adjustl(scalingtype)))
  if (.not. isok) then
     write(*, *) 'Error parsing parameter scaling_type, stopping!'
     stop
  end if
  call splitstr(scalingtype, ',', scalingtypes)
end subroutine parsescalingtype

!> Writes the parameters and observables of this module and its
!! dependencies. The @p writer defines the format and output unit..
subroutine mc_sweep_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  character(len=200) :: joined
  call writecomment(writer, 'MC sweep parameters')
  call writeparameter(writer, 'move_ratio', moveratio)
  call writeparameter(writer, 'scaling_ratio', scalingratio)
  call writeparameter(writer, 'max_scaling', maxscaling)
  call writeparameter(writer, 'pt_period', pt_period)
  if (allocated(scalingtypes)) then
     call join(scalingtypes, ',', joined)
  else
     joined = ''
  end if
  call writeparameter(writer, 'scaling_type', trim(joined))
  call writeparameter(writer, 'nmovetrials', nmovetrials)
  call writeparameter(writer, 'nacceptedmoves', nacceptedmoves)
  call writeparameter(writer, 'current_move_ratio', &
       real(nacceptedmoves, dp)/real(nmovetrials, dp))
  call writeparameter(writer, 'nscalingtrials', nscalingtrials)
  call writeparameter(writer, 'nacceptedscalings', nacceptedscalings)
  call writeparameter(writer, 'current_scaling_ratio', &
       real(nacceptedscalings, dp)/real(nscalingtrials, dp))
  call writeparameter(writer, 'pressure', pressure)
  call writeparameter(writer, 'temperature', temperature)
  call writeparameter(writer, 'volume', currentvolume)
  call writeparameter(writer, 'enthalpy', etotal + currentvolume * pressure)
  call writeparameter(writer, 'total_energy', etotal)
  call particlemover_writeparameters(writer)
  call energy_writeparameters(writer)
end subroutine mc_sweep_writeparameters

!> Setter for @p the_particles and the simulation box @p the_simbox.
subroutine set_system(the_simbox, the_particles)
  type(poly_box), intent(in) :: the_simbox
  type(particledat), intent(in) :: the_particles(:)
  logical :: overlap
  logical :: should_allocate
  !real(dp) :: debug_etotal
  !logical :: debug_overlap 
  should_allocate = .false.
  simbox = the_simbox
  if(allocated(particles)) deallocate(particles)
  allocate(particles(size(the_particles)))
  particles = the_particles
  call update(sl, simbox, particles)
  call totalenergy(sl, simbox, particles, etotal, overlap)
  if (overlap) stop 'mc_sweep:set_system: Trying to set a geometry with overlap!' 
  currentvolume = volume(simbox)
end subroutine set_system

!> Getter for @p the_particles and the simulation box @p the_simbox.
subroutine get_system(the_simbox, the_particles)
  type(poly_box), intent(out) :: the_simbox
  type(particledat), allocatable, intent(out) :: the_particles(:)
  the_simbox = simbox
  allocate(the_particles(size(particles)))
  the_particles = particles
end subroutine get_system

!> Runs one sweep of Metropolis Monte Carlo updates to the system. A
!! full Parallel tempering NPT-ensemble sweep consists of trial moves
!! of particles, trial scaling of the simulation box (barostat) and an
!! exchange of particle system coordinates with another particle system
!! (replica) in another temperature (replica exchange).
!! 
!! @param genstates random number generator states for all threads.
!! @param isweep the sweep counter.
!!  
subroutine sweep(genstates, isweep)    
  type(rngstate), intent(inout) :: genstates(0:)
  integer, intent(in) :: isweep
  integer :: ivolmove
  real(dp) :: beta
  call update(sl, simbox, particles)
  call make_particle_moves(simbox, particles, genstates, sl)
  call update(sl, simbox, particles)
  do ivolmove = 1, size(scalingtypes)
     call movevol(simbox, particles, scalingtypes(ivolmove), genstates(0))
  end do
  if (mod(isweep, pt_period) == 0) then
     beta = 1._dp/temperature
     call try_beta_exchanges(beta, etotal, 3, genstates(0)) 
     temperature = 1._dp/beta
  end if
end subroutine sweep

!> Schedules parallel moves of particles using OpenMP with a domain 
!! decomposition algorithm. 
!!
!! @see e.g. G. Heffelfinger and M. Lewitt. J. Comp. Chem., 17(2):250–265,
!! 1996.
!! 
!! @param simbox is the simulation box in which the particles reside.
!! @param particles is the array of particles to move.
!! @param genstates random number generator states. Each thread needs one.
!! @param sl the cell list presenting the decomposition.
!! 
subroutine make_particle_moves(simbox, particles, genstates, sl)
  implicit none
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(inout) :: particles(:)
  type(rngstate), intent(inout) :: genstates(0:)
  type(simplelist), intent(in) :: sl
  !$ integer :: n_threads = 1
  integer :: thread_id = 0
  real(dp) :: dE 
  logical :: isaccepted
  !! You may be tempted to make the allocatable arrays automatic, but there's
  !! no performance gained and depending on the system large automatic arrays
  !! may cause a stack overflow.
  type(particledat), allocatable :: temp_particles(:)
  logical, allocatable :: nbr_mask(:)
  integer :: n_cell ! particles in the cell where particles are moved.
  integer :: n_local ! particles cell and its neighbour cells.
  real(dp) :: dE_ij
  integer :: nacc, ntrials
  integer :: j, ix, iy, iz, jx, jy, jz
  dE = 0._dp
  dE_ij = 0._dp
  nacc = 0
  ntrials = 0
  !! Loop over cells. This can be thought of as looping through a 
  !! 2 x 2 x 2 cube of cells.
  !$OMP PARALLEL shared(particles, simbox, sl, genstates)& 
  !$OMP& private(thread_id, n_threads, temp_particles, nbr_mask)&
  !$OMP& reduction(+:dE, nacc, ntrials) 
  !$ thread_id = omp_get_thread_num()
  allocate(temp_particles(minval([size(particles),&
       27 * maxval(sl%counts)])), nbr_mask(size(particles)))
  do iz=0, min(1, sl%nz-1)
     do iy=0, min(1, sl%ny-1)
        do ix=0, min(1, sl%nx-1)
           !$OMP DO collapse(3) private(j, n_cell, n_local, dE_ij, isaccepted)&
           !$OMP& schedule(dynamic)
           do jz = iz, sl%nz - 1, 2
              do jy = iy, sl%ny - 1, 2
                 do jx = ix, sl%nx - 1, 2
                    nbr_mask = .false.
                    n_cell = sl%counts(jx, jy, jz)
                    call simplelist_nbrmask(sl, simbox, jx, jy, jz, nbr_mask)
                    n_local = count(nbr_mask)
                    nbr_mask(sl%indices(1:n_cell, jx, jy, jz)) = .false.
                    temp_particles(1:n_cell) = &
                         particles(sl%indices(1:n_cell, jx, jy, jz))
                    temp_particles(n_cell + 1 : n_local) = &
                         pack(particles, nbr_mask)
                    do j = 1, n_cell
                       call moveparticle_2(simbox, &
                            temp_particles(1:n_local), j, &
                        genstates(thread_id), dE_ij, isaccepted)
                       if (isaccepted) then 
                          nacc = nacc + 1
                          dE = dE + dE_ij
                       end if
                       ntrials = ntrials + 1
                    end do
                    particles(sl%indices(1:n_cell, jx, jy, jz)) = &
                         temp_particles(1:n_cell)
                 end do
              end do
           end do
           !$OMP END DO 
           !! The end of parallelized loop forces an implicit barrier.
           !! Memory view of the threads is also synchronized here.
        end do
        !$OMP BARRIER
     end do
     !$OMP BARRIER
  end do
  if (allocated(temp_particles)) deallocate(temp_particles)
  if (allocated(nbr_mask)) deallocate(nbr_mask)
  !$OMP END PARALLEL
  etotal = etotal + dE
  nmovetrials = nmovetrials + ntrials
  nacceptedmoves = nacceptedmoves + nacc
end subroutine make_particle_moves

!> Performs a trial move of @p particles(@p i). @p simbox is the
!! simulation box in which the @p particles reside. @p genstate is the
!! random number generator state. After the move, @p dE contains the
!! change in energy of the system and @p isaccepted == .true. if the
!! move was accepted. 
pure subroutine moveparticle_2(simbox, particles, i, genstate, dE, isaccepted)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(inout) :: particles
  integer, intent(in) :: i
  type(rngstate), intent(inout) :: genstate
  real(dp), intent(out) :: dE
  logical, intent(out) :: isaccepted
  
  type(particledat) :: newparticle
  type(particledat) :: oldparticle
  logical :: overlap
  real(dp) :: enew
  real(dp) :: eold
  
  enew = 0._dp
  eold = 0._dp
  dE = 0._dp
  overlap = .false.
  isaccepted = .false.
  newparticle = particles(i)
  call move(newparticle, genstate)
  call setposition(newparticle, minimage(simbox, position(newparticle)))
  oldparticle = particles(i) 
  particles(i) = newparticle
  call singleparticleenergy(simbox, particles, i, enew, overlap)
  
  particles(i) = oldparticle
  if(.not. overlap) then 
     call singleparticleenergy(simbox, particles, i, eold, overlap)
     call acceptchange(eold, enew, genstate, isaccepted)
     if(isaccepted) then
        particles(i) = newparticle
        dE = enew - eold
     end if
  end if  
end subroutine moveparticle_2

!> Returns the simulation temperature.
pure function gettemperature() result(temp)
  real(dp) :: temp
  temp = temperature
end function gettemperature

!> Sets the simulation temperature.
subroutine settemperature(temperaturein)
  real(dp), intent(in) :: temperaturein
  temperature = temperaturein
end subroutine settemperature

!> Performs a trial volume scaling which scales the @p simbox and all
!! the positions of @p particles. For more information see for example
!! Allen and Tildesley: Computer Simulation of Liquids: the chapter
!! about NPT ensemble Monte Carlo.
!! 
!! @param simbox the simulation box.
!! @param particles the particles in the simulation box.
!! @param scalingtype string defining the direction(s) for changing
!!        the system volume.
!! @param genstate the random number generator state.
!! 
subroutine movevol(simbox, particles, scalingtype, genstate)    
  type(poly_box), intent(inout) :: simbox
  type(particledat), dimension(:), intent(inout) :: particles
  character(len=*), intent(in) :: scalingtype
  type(rngstate), intent(inout) :: genstate
  integer :: nparticles 
  logical :: overlap
  real(dp) :: Vo, Vn
  logical :: isaccepted
  real(dp) :: totenew
  real(dp) :: boltzmannn
  real(dp) :: boltzmanno
  type(poly_box) :: oldbox
  real(dp), dimension(3) :: scaling
  logical :: is_update_needed
  nparticles = size(particles)
  overlap = .false.
  is_update_needed = .false.
  
  !! Store old volume and simulation box
  Vo = volume(simbox)
  oldbox = simbox
  
  !! It seems that total energy may drift if it is not updated here:
  call totalenergy(sl, simbox, particles, etotal, overlap)
  if (overlap) stop 'movevol: overlap in old configuration! '//&
       'Should never happen!'
  
  !! Scale coordinates and the simulations box
  scaling = genvoltrial_scale(simbox, maxscaling, genstate, &
       trim(adjustl(scalingtype)))
  call scalepositions(oldbox, simbox, particles, nparticles) 
  Vn = volume(simbox)
  
  !! Check that new dimensions are ok. 
  call check_simbox(simbox)

  if (Vn/Vo < 1._dp/(1._dp + 2._dp * get_max_translation()/get_cutoff())) then
     !! Volume shrank so quickly that the neighbourlist needs updating.
     call update(sl, simbox, particles)
     is_update_needed = .true.
  end if
  
    !! Calculate potential energy in the scaled system.
  call totalenergy(sl, simbox, particles, totenew, overlap)
  
  if (overlap) then
     !! Scale particles back to old coordinates.
     call scalepositions(simbox, oldbox, particles, nparticles)
     simbox = oldbox
     if (is_update_needed) call update(sl, simbox, particles)
  else
     boltzmannn = totenew + pressure * Vn - real(nparticles, dp) * &
          temperature * log(Vn)  
     boltzmanno = etotal + pressure * Vo - real(nparticles, dp) * &
          temperature * log(Vo)
     call acceptchange(boltzmanno, boltzmannn, genstate, isaccepted)
     if (isaccepted) then
        etotal = totenew
        nacceptedscalings = nacceptedscalings + 1
        currentvolume = volume(simbox)
        call update(sl, simbox, particles)
     else 
        !! Scale particles back to old coordinates
        call scalepositions(simbox, oldbox, particles, nparticles)
        simbox = oldbox
        if (is_update_needed) call update(sl, simbox, particles)
     end if
  end if
  nscalingtrials = nscalingtrials + 1
end subroutine movevol

!> Check that @p simbox is large enough if it is periodic.
subroutine check_simbox(simbox)
  type(poly_box), intent(in) :: simbox
  if (simbox%xperiodic .and. getx(simbox) < 2._dp * get_cutoff()) &
       stop 'Simulation box too small!'
  if (simbox%yperiodic .and. gety(simbox) < 2._dp * get_cutoff()) &
       stop 'Simulation box too small!'
  if (simbox%zperiodic .and. getz(simbox) < 2._dp * get_cutoff()) &
       stop 'Simulation box too small!'
end subroutine

!> Implements the Metropolis acceptance rule for a Monte Carlo update. 
!!
!! @param oldenergy the energy/enthalpy of the system before the move.
!! @param newenergy the energy/enthalpy of the system after the move.
!! @param genstate is the random number generator state.
!! @param isaccepted == .true. if the move is accepted. 
!!
pure subroutine acceptchange(oldenergy, newenergy, genstate, isaccepted)
  real(dp), intent(in) :: oldenergy
  real(dp), intent(in) :: newenergy
  type(rngstate), intent(inout) :: genstate
  logical, intent(out) :: isaccepted
  real(dp) :: dE
  real(dp) :: r
  isaccepted = .true.
  dE = newenergy - oldenergy
  if (dE > 0._dp) then
     call rng(genstate, r)
     isaccepted = (r < exp(-dE/temperature))
  end if
end subroutine acceptchange

!> Adjusts the maximum values for trial moves of particles and trial
!! scalings of the simulation volume. Should be used only during
!! equilibration run. 
subroutine updatemaxvalues
  real(dp) :: newdthetamax
  real(dp) :: newdximax
  !! Adjust scaling
  maxscaling = newmaxvalue(nscalingtrials, nacceptedscalings, scalingratio,&
       maxscaling)
  
  call getmaxmoves(newdximax, newdthetamax)
  !! Adjust translation
  newdximax = newmaxvalue(nmovetrials, nacceptedmoves, moveratio, newdximax)
  
  !! Update the minimum cell side length of the cell list because the maximum
  !! translation has changed: 
  sl%min_length = get_cutoff() + 2._dp * get_max_translation()
  
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
  
!> Returns the simulation pressure in reduced units.
function getpressure()
  real(dp) :: getpressure
  getpressure = pressure
end function getpressure

!> Resets the counters that are used to monitor acceptances of trial
!! moves. 
subroutine resetcounters
  nacceptedmoves = 0
  nmovetrials = 0
  nscalingtrials = 0 
  nacceptedscalings = 0
end subroutine resetcounters

!> Scales the positions of @p particles with the same factors that were
!! used to scale the simulation box dimensions from @p oldbox to
!! @p newbox.
subroutine scalepositions(oldbox, newbox, particles, nparticles)
  type(poly_box), intent(in) :: oldbox
  type(poly_box), intent(in) :: newbox
  type(particledat), dimension(:), intent(inout) :: particles
  integer, intent(in) :: nparticles
  integer :: i
  do i = 1, nparticles
     particles(i)%x = particles(i)%x * getx(newbox) / getx(oldbox)
     particles(i)%y = particles(i)%y * gety(newbox) / gety(oldbox)
     particles(i)%z = particles(i)%z * getz(newbox) / getz(oldbox)
  end do
end subroutine scalepositions

!> Tests that there are no overlaps in the current configuration of
!! particles and the simulation box.
subroutine test_configuration()
  real(dp) :: total_e
  logical :: overlap
  if (.not. is_initialized) then
     stop 'mc_sweep:test_configuration: Error: module not initialized!'
  end if
  !! This is pretty heavy since goes through all particles:
  call totalenergy(simbox, particles, total_e, overlap)
  if (overlap) then
     stop 'mc_sweep:test_configuration: Overlap!'
  end if
end subroutine test_configuration


!include 'map_and_reduce.f90'

end module mc_sweep
