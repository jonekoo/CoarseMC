!> Implements Metropolis Monte Carlo updates of the coordinates of the
!! particles and the simulation box dimensions. The module also
!! controls the parallel tempering updates. If compiled and run
!! with OpenMP, the domain decomposition algorithm is used for the
!! particle moves and energy calculations. 
module mc_sweep
  use nrtype
  use class_poly_box
  use particle
  use class_parameterizer
  use class_parameter_writer
  use genvoltrial
  use utils 
  !$ use omp_lib
  use class_simplelist
  include 'rng.inc'
  use m_particlegroup, only: particlegroup
  implicit none  

  !> The maximum absolute change of volume in a trial volume update.
  real(dp), save :: maxscaling = 100._dp

  !> The total energy.
  real(dp), save :: etotal = 0._dp

  !> The types of trial volume updates.
  character(len = 200), dimension(:), allocatable, save :: scalingtypes

  !> True if the module is correctly initialized.
  logical :: is_initialized = .false.

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
subroutine mcsweep_init(reader)
  type(parameterizer), intent(in) :: reader
  character(len = 200), save :: scalingtype = "z"
  call getparameter(reader, 'scaling_type', scalingtype)
  call parsescalingtype(scalingtype)
  call getparameter(reader, 'max_scaling', maxscaling)
  is_initialized = .true.
end subroutine mcsweep_init

!> Finalizes the module state.
subroutine mcsweep_finalize
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

subroutine set_maxscaling(upper_limit)
  real(dp), intent(in) :: upper_limit
  maxscaling = upper_limit
end subroutine set_maxscaling

function get_maxscaling() result(upper_limit)
  real(dp) :: upper_limit
  upper_limit = maxscaling
end function get_maxscaling

!> Writes the parameters and observables of this module and its
!! dependencies. The @p writer defines the format and output unit..
subroutine mc_sweep_writeparameters(writer)
  type(parameter_writer), intent(inout) :: writer
  character(len=200) :: joined
  call writecomment(writer, 'MC sweep parameters')
  call writeparameter(writer, 'max_scaling', maxscaling)
  if (allocated(scalingtypes)) then
     call join(scalingtypes, ',', joined)
  else
     joined = ''
  end if
  call writeparameter(writer, 'scaling_type', trim(joined))
  call writeparameter(writer, 'total_energy', etotal)
end subroutine mc_sweep_writeparameters

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
subroutine make_particle_moves(simbox, group, genstates, subr_particle_energy, &
     temperature, n_trials, n_accepted)
  implicit none
  type(poly_box), intent(in) :: simbox
  type(particlegroup), intent(inout) :: group
  type(rngstate), intent(inout) :: genstates(0:)
  procedure(particle_energy) :: subr_particle_energy
  real(dp), intent(in) :: temperature
  integer, intent(out) :: n_accepted, n_trials
  
  !$ integer :: n_threads
  integer :: thread_id
  real(dp) :: dE 
  logical :: isaccepted
  !! You may be tempted to make the allocatable arrays automatic, but there's
  !! no performance gained and depending on the system large automatic arrays
  !! may cause a stack overflow.
  type(particledat), allocatable :: temp_particles(:)
  logical, allocatable :: nbr_mask(:)
  integer :: n_cell ! particles in the cell where particles are moved.
  integer :: n_local ! particles cell and its neighbour cells.
  real(dp) :: dE_j
  integer :: j, ix, iy, iz, jx, jy, jz, n_trials_j, n_accepted_j
  call update(group%sl, simbox, group%particles)
  thread_id = 0
  !$ n_threads = 1
  dE = 0._dp
  dE_j = 0._dp
  n_accepted = 0
  n_trials = 0
  !! Loop over cells. This can be thought of as looping through a 
  !! 2 x 2 x 2 cube of cells.
  !$OMP PARALLEL shared(group, simbox, genstates)& 
  !$OMP& private(thread_id, n_threads, temp_particles, nbr_mask, n_trials_j, &
  !$OMP& n_accepted_j)&
  !$OMP& reduction(+:dE, n_accepted, n_trials) 
  !$ thread_id = omp_get_thread_num()
  allocate(temp_particles(minval([size(group%particles),&
       27 * maxval(group%sl%counts)])), nbr_mask(size(group%particles)))
  do iz=0, min(1, group%sl%nz-1)
     do iy=0, min(1, group%sl%ny-1)
        do ix=0, min(1, group%sl%nx-1)
           !$OMP DO collapse(3) private(j, n_cell, n_local, dE_j, isaccepted)&
           !$OMP& schedule(dynamic)
           do jz = iz, group%sl%nz - 1, 2
              do jy = iy, group%sl%ny - 1, 2
                 do jx = ix, group%sl%nx - 1, 2
                    nbr_mask = .false.
                    n_cell = group%sl%counts(jx, jy, jz)
                    call simplelist_nbrmask(group%sl, simbox, jx, jy, jz, &
                         nbr_mask)
                    n_local = count(nbr_mask)
                    nbr_mask(group%sl%indices(1:n_cell, jx, jy, jz)) = .false.
                    temp_particles(1:n_cell) = &
                         group%particles(group%sl%indices(1:n_cell, jx, jy, jz))
                    temp_particles(n_cell + 1 : n_local) = &
                         pack(group%particles, nbr_mask)
                    do j = 1, n_cell
                       call moveparticle_2(simbox, &
                            temp_particles(1:n_local), j, temperature, &
                            genstates(thread_id), subr_particle_energy, dE_j, &
                            n_trials=n_trials_j, n_accepted=n_accepted_j)
                       n_accepted = n_accepted + n_accepted_j
                       n_trials = n_trials + n_trials_j
                       dE = dE + dE_j
                    end do
                    group%particles(group%sl%indices(1:n_cell, jx, jy, jz)) = &
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
end subroutine make_particle_moves

!> Returns the simulation temperature.
!pure function gettemperature() result(temp)
!  real(dp) :: temp
!  temp = temperature
!end function gettemperature

!> Sets the simulation temperature.
!subroutine settemperature(temperaturein)
!  real(dp), intent(in) :: temperaturein
!  temperature = temperaturein
!end subroutine settemperature


subroutine update_volume(simbox, group, genstate, temperature, pressure, &
     subr_particle_energy, n_trials, n_accepted)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup), intent(inout) :: group
  type(rngstate), intent(inout) :: genstate
  real(dp), intent(in) :: temperature, pressure
  procedure(particle_energy) :: subr_particle_energy
  integer, intent(out), optional :: n_trials, n_accepted
  integer :: i
  do i = 1, size(scalingtypes)
     call movevol(simbox, group, scalingtypes(i), genstate, temperature, &
          pressure, subr_particle_energy, n_trials, n_accepted)
  end do
end subroutine update_volume


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
subroutine movevol(simbox, group, scalingtype, genstate, temperature, pressure,&
     subr_particle_energy, n_trials, n_accepted)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup), intent(inout) :: group
  character(len=*), intent(in) :: scalingtype
  type(rngstate), intent(inout) :: genstate
  real(dp), intent(in) :: temperature, pressure
  procedure(particle_energy) :: subr_particle_energy
  integer, intent(out), optional :: n_trials, n_accepted
  integer :: nparticles 
  logical :: overlap
  real(dp) :: Vo, Vn
  logical :: isaccepted
  real(dp) :: totenew
  real(dp) :: boltzmannn
  real(dp) :: boltzmanno
  type(poly_box) :: oldbox
  real(dp), dimension(3) :: scaling
  nparticles = size(group%particles)
  overlap = .false.
  
  !! Store old volume and simulation box
  Vo = volume(simbox)
  oldbox = simbox
  
  !! It seems that total energy may drift if it is not updated here:
  call total_energy(group%sl, simbox, group%particles, subr_particle_energy, &
       etotal, overlap)
  if (overlap) stop 'movevol: overlap in old configuration! '//&
       'Should never happen!'
  
  !! Scale coordinates and the simulations box
  scaling = genvoltrial_scale(simbox, maxscaling, genstate, &
       trim(adjustl(scalingtype)))
  call scalepositions(oldbox, simbox, group%particles, nparticles) 
  Vn = volume(simbox)
  
  !! Calculate potential energy in the scaled system.
  call total_energy(group%sl, simbox, group%particles, subr_particle_energy, &
       totenew, overlap)
  
  if (overlap) then
     !! Scale particles back to old coordinates.
     call scalepositions(simbox, oldbox, group%particles, nparticles)
     simbox = oldbox
  else
     boltzmannn = totenew + pressure * Vn - real(nparticles, dp) * &
          temperature * log(Vn)  
     boltzmanno = etotal + pressure * Vo - real(nparticles, dp) * &
          temperature * log(Vo)
     call acceptchange(boltzmanno, boltzmannn, temperature, genstate, &
          isaccepted)
     if (isaccepted) then
        etotal = totenew
     else 
        !! Scale particles back to old coordinates
        call scalepositions(simbox, oldbox, group%particles, nparticles)
        simbox = oldbox
     end if
  end if

  if (present(n_accepted)) then
     if (isaccepted) then
        n_accepted = 1
     else
        n_accepted = 0
     end if
  end if

  if(present(n_trials)) n_trials = 1
end subroutine movevol

!> Computes the total @p energy of the system. Interactions are computed
!! cell-by-cell for the @p particles in the cell list @p sl. @p simbox
!! is the simulation cell. If any two particles are too close to each
!! other, @p overlap is true. 
subroutine total_energy(sl, simbox, particles, subr_particle_energy, &
     energy, overlap)
  type(simplelist), intent(inout) :: sl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  procedure(particle_energy) :: subr_particle_energy
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap 
  integer :: i, j, ix, iy, iz
  logical :: mask(size(particles))
  real(dp) :: energy_j
  logical :: overlap_j
  integer :: helper(size(particles))
  integer, allocatable :: temp_helper(:)
  integer :: n_mask
  integer :: temp_j
  type(particledat), allocatable :: temp_particles(:)

  call update(sl, simbox, particles)
  helper = (/(i, i=1, size(particles))/) !! ifort vectorizes
  energy = 0._dp
  overlap = .false.
  
  !$OMP PARALLEL default(shared) reduction(+:energy) reduction(.or.:overlap)& 
  !$OMP& private(energy_j, overlap_j, i, j, mask, temp_particles, temp_j, temp_helper, n_mask)
  allocate(temp_particles(size(particles)), temp_helper(size(particles))) 
  !$OMP DO collapse(3) schedule(dynamic)
  do ix=0, sl%nx-1 
  do iy=0, sl%ny-1
  do iz=0, sl%nz-1
    call simplelist_nbrmask(sl, simbox, ix, iy, iz, mask)
    n_mask = count(mask) 
    temp_particles(1:n_mask) = pack(particles, mask)
    temp_helper(1:n_mask) = pack(helper, mask)
    do i=1, sl%counts(ix, iy, iz)
      j = sl%indices(i, ix, iy, iz)
      ! Find position of particles(j) in temp_particles:
      do temp_j = 1, n_mask 
         if(temp_helper(temp_j) == j) exit
      end do
      call subr_particle_energy(simbox, temp_particles(temp_j:n_mask), &
           1, energy_j, overlap_j)
      overlap = overlap .or. overlap_j 
      energy = energy + energy_j
    end do 
  end do
  end do
  end do
  !$OMP END DO 
  if (allocated(temp_particles)) deallocate(temp_particles)
  if (allocated(temp_helper)) deallocate(temp_helper)
  !$OMP END PARALLEL  
end subroutine total_energy


!> Returns the simulation pressure in reduced units.
!function getpressure()
!  real(dp) :: getpressure
!  getpressure = pressure
!end function getpressure

!> Scales the positions of @p particles with the same factors that were
!! used to scale the simulation box dimensions from @p oldbox to
!! @p newbox.
pure subroutine scalepositions(oldbox, newbox, particles, nparticles)
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
!subroutine test_configuration(simbox, particles)
!  real(dp) :: total_e
!  logical :: overlap
!  type(poly_box), intent(in) :: simbox
!  type(particledat), intent(in) :: particles(:)
!  if (.not. is_initialized) then
!     stop 'mc_sweep:test_configuration: Error: module not initialized!'
!  end if
!  !! This is pretty heavy since goes through all particles:
!  call total_energy(simbox, particles, total_e, overlap)
!  if (overlap) then
!     stop 'mc_sweep:test_configuration: Overlap!'
!  end if
!end subroutine test_configuration


!include 'map_and_reduce.f90'

end module mc_sweep
