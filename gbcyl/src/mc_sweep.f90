!> Implements Metropolis Monte Carlo updates of the coordinates of the
!! particles and the simulation box dimensions. The module also
!! controls the parallel tempering updates. If compiled and run
!! with OpenMP, the domain decomposition algorithm is used for the
!! particle moves and energy calculations. 
module mc_sweep
  use nrtype
  use energy, only: get_cutoff, energy_init, energy_writeparameters
  use class_poly_box
  use particle
  use particle_mover, only: get_max_translation, getmaxmoves, &
       particlemover_init, particlemover_writeparameters, setmaxmoves
  use class_parameterizer
  use class_parameter_writer
  use genvoltrial
  use utils 
  !$ use omp_lib
  use class_simplelist
  include 'rng.inc'
  use m_particlegroup, only: particlegroup
  implicit none  

  !> The number of accepted particle moves
  integer, save :: nacceptedmoves

  !> The number of accepted volume moves
  integer, save :: nacceptedscalings

  !> The maximum absolute change of volume in a trial volume update.
  real(dp), save :: maxscaling = 100._dp

  !> The desired acceptance ratio for trial particle moves.
  real(dp), save :: moveratio

  !> The desired acceptance ratio for trial volume updates.
  real(dp), save :: scalingratio

  !> The total energy.
  real(dp), save :: etotal = 0._dp

  !> Current volume of the system.
  !real(dp), save :: currentvolume = 0._dp

  !> Counter for trial particle moves.
  integer, save :: nmovetrials = 0

  !> Counter for trial volume updates.
  integer, save :: nscalingtrials = 0

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
  !type(particlegroup), intent(inout) :: group
  real(dp) :: min_cell_length
  character(len = 200), save :: scalingtype = "z"
  call particlemover_init(reader)
  call energy_init(reader)
  call getparameter(reader, 'scaling_type', scalingtype)
  call parsescalingtype(scalingtype)
  call getparameter(reader, 'move_ratio', moveratio)
  call getparameter(reader, 'scaling_ratio', scalingratio)
  call getparameter(reader, 'max_scaling', maxscaling)
  call getparameter(reader, 'nmovetrials', nmovetrials)
  call getparameter(reader, 'nscalingtrials', nscalingtrials)
  call getparameter(reader, 'nacceptedscalings', nacceptedscalings)
  !! The initializations below may affect restart so that it does not result
  !! in the same simulation. 
  nacceptedmoves = 0
  nacceptedscalings = 0
  !min_cell_length = get_cutoff() + 2._dp * get_max_translation()
  !!call simplelist_init(reader)
  !!$ if (.true.) then
  !!$ call new_simplelist(sl, simbox, particles, min_cell_length, &
  !!$& is_x_even = isxperiodic(simbox), is_y_even = isyperiodic(simbox), &
  !!$& is_z_even = iszperiodic(simbox), cutoff=get_cutoff())
  !!$ else 
  !call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
  !     cutoff=get_cutoff())
  !!$ end if
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

!> Writes the parameters and observables of this module and its
!! dependencies. The @p writer defines the format and output unit..
subroutine mc_sweep_writeparameters(writer)
  type(parameter_writer), intent(inout) :: writer
  character(len=200) :: joined
  call writecomment(writer, 'MC sweep parameters')
  call writeparameter(writer, 'move_ratio', moveratio)
  call writeparameter(writer, 'scaling_ratio', scalingratio)
  call writeparameter(writer, 'max_scaling', maxscaling)
  if (allocated(scalingtypes)) then
     call join(scalingtypes, ',', joined)
  else
     joined = ''
  end if
  call writeparameter(writer, 'scaling_type', trim(joined))
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
  call writeparameter(writer, 'total_energy', etotal)
  call particlemover_writeparameters(writer)
  call energy_writeparameters(writer)
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
     temperature)
  implicit none
  type(poly_box), intent(in) :: simbox
  type(particlegroup), intent(inout) :: group
  type(rngstate), intent(inout) :: genstates(0:)
  procedure(particle_energy) :: subr_particle_energy
  real(dp), intent(in) :: temperature
  
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
  real(dp) :: dE_ij
  integer :: nacc, ntrials
  integer :: j, ix, iy, iz, jx, jy, jz
  call update(group%sl, simbox, group%particles)
  thread_id = 0
  !$ n_threads = 1
  dE = 0._dp
  dE_ij = 0._dp
  nacc = 0
  ntrials = 0
  !! Loop over cells. This can be thought of as looping through a 
  !! 2 x 2 x 2 cube of cells.
  !$OMP PARALLEL shared(group, simbox, genstates)& 
  !$OMP& private(thread_id, n_threads, temp_particles, nbr_mask)&
  !$OMP& reduction(+:dE, nacc, ntrials) 
  !$ thread_id = omp_get_thread_num()
  allocate(temp_particles(minval([size(group%particles),&
       27 * maxval(group%sl%counts)])), nbr_mask(size(group%particles)))
  do iz=0, min(1, group%sl%nz-1)
     do iy=0, min(1, group%sl%ny-1)
        do ix=0, min(1, group%sl%nx-1)
           !$OMP DO collapse(3) private(j, n_cell, n_local, dE_ij, isaccepted)&
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
                            genstates(thread_id), subr_particle_energy, dE_ij, &
                            isaccepted)
                       if (isaccepted) then 
                          nacc = nacc + 1
                          dE = dE + dE_ij
                       end if
                       ntrials = ntrials + 1
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
  nmovetrials = nmovetrials + ntrials
  nacceptedmoves = nacceptedmoves + nacc
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
     subr_particle_energy)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup), intent(inout) :: group
  type(rngstate), intent(inout) :: genstate
  real(dp), intent(in) :: temperature, pressure
  procedure(particle_energy) :: subr_particle_energy
  integer :: i
  do i = 1, size(scalingtypes)
     call movevol(simbox, group, scalingtypes(i), genstate, temperature, &
          pressure, subr_particle_energy)
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
     subr_particle_energy)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup), intent(inout) :: group
  character(len=*), intent(in) :: scalingtype
  type(rngstate), intent(inout) :: genstate
  real(dp), intent(in) :: temperature, pressure
  procedure(particle_energy) :: subr_particle_energy
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
  nparticles = size(group%particles)
  overlap = .false.
  call update(group%sl, simbox, group%particles)
  is_update_needed = .false.
  
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
  
  !! Check that new dimensions are ok. 
  call check_simbox(simbox)

  if (Vn/Vo < 1._dp/(1._dp + 2._dp * get_max_translation()/get_cutoff())) then
     !! Volume shrank so quickly that the neighbourlist needs updating.
     call update(group%sl, simbox, group%particles)
     is_update_needed = .true.
  end if
  
    !! Calculate potential energy in the scaled system.
  call total_energy(group%sl, simbox, group%particles, subr_particle_energy, &
       totenew, overlap)
  
  if (overlap) then
     !! Scale particles back to old coordinates.
     call scalepositions(simbox, oldbox, group%particles, nparticles)
     simbox = oldbox
     if (is_update_needed) call update(group%sl, simbox, group%particles)
  else
     boltzmannn = totenew + pressure * Vn - real(nparticles, dp) * &
          temperature * log(Vn)  
     boltzmanno = etotal + pressure * Vo - real(nparticles, dp) * &
          temperature * log(Vo)
     call acceptchange(boltzmanno, boltzmannn, temperature, genstate, &
          isaccepted)
     if (isaccepted) then
        etotal = totenew
        nacceptedscalings = nacceptedscalings + 1
        call update(group%sl, simbox, group%particles)
     else 
        !! Scale particles back to old coordinates
        call scalepositions(simbox, oldbox, group%particles, nparticles)
        simbox = oldbox
        if (is_update_needed) call update(group%sl, simbox, group%particles)
     end if
  end if
  nscalingtrials = nscalingtrials + 1
end subroutine movevol

!> Computes the total @p energy of the system. Interactions are computed
!! cell-by-cell for the @p particles in the cell list @p sl. @p simbox
!! is the simulation cell. If any two particles are too close to each
!! other, @p overlap is true. 
subroutine total_energy(sl, simbox, particles, subr_particle_energy, &
     energy, overlap)
  type(simplelist), intent(in) :: sl
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

!> Adjusts the maximum values for trial moves of particles and trial
!! scalings of the simulation volume. Should be used only during
!! equilibration run. 
subroutine updatemaxvalues(group)
  type(particlegroup), intent(inout) :: group
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
!function getpressure()
!  real(dp) :: getpressure
!  getpressure = pressure
!end function getpressure

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
