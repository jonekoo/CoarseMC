module mc_sweep
  use nrtype
  use energy
  use class_poly_box
  use particle
  use class_poly_nbrlist
  use pt
  use class_parameterizer
  use class_parameter_writer
  use mpi
  use genvoltrial
  use utils 
  use gayberne
  !$ use omp_lib
  !$ use class_simplelist
  include 'rng.inc'
  implicit none  
  private 

  public :: init
  public :: sweep
  public :: updatemaxvalues
  public :: getpressure
  public :: gettemperature
  public :: mc_sweep_writeparameters
  public :: settemperature
  public :: resetcounters
  public :: movevol
  public :: makeptmove
  public :: moveparticle
  public :: volratio
  public :: ptratio
  public :: getscalingtypes

  integer, save :: nacceptedmoves
  !! number of accepted particle moves

  integer, save :: nacceptedscalings
  !! number of accepted volume moves

  real(dp), save :: maxscaling
  !! the maximum absolute scaling that can happen in a volume move
  !! currently redundant parameter for selecting the type of volume scaling

  real(dp), save :: temperature
  !! the simulation temperature that is used in the Metropolis acceptance
  !! criterion

  real(dp), save :: pressure
  !! the simulation pressure that is used in the Metropolis acceptance
  !! criterion for trial volume scalings.

  integer, save :: adjusttype
  !! currently redundant parameter for selecting the way of adjusting the 
  !! maximum translation and maximum rotation parameters

  real(dp), save :: moveratio
  !! the desired acceptance ratio for trial particle moves

  real(dp), save :: scalingratio
  !! the desired acceptance ratio for trial volume scalings

  real(dp), save :: largesttranslation = 1._dp    
  !! largest translation that can happen. Default is 1 x molecule diameter 

  real(dp), save :: smallesttranslation = 0.1_dp 
  !! The maxtranslation won't be adjusted below this. Default is 1/10 molecule
  !! diameter.

  real(dp), save :: largestrotation = 1.57_dp     !! ~ pi/4
  real(dp), save :: smallestrotation = 0.157_dp  !! 1/10 largest rotation
  real(dp), save :: ptlow
  real(dp), save :: pthigh
  integer, save :: ptratio = 1
  integer, save :: volratio = -1
  integer, save :: radialratio = -1
  real(dp), save :: etotal = 0._dp
  real(dp), save :: currentvolume = 0._dp
  integer, save :: nmovetrials = 0
  integer, save :: nscalingtrials = 0
  integer, save :: nradialtrials = 0
  integer, save :: nacceptedradial = 0
  type(poly_nbrlist), pointer, save :: nbrlist
  character(len = 200), save :: scalingtype = "z"
  character(len = 200), dimension(:), pointer, save :: scalingtypes

  interface init
    module procedure initparameterizer
  end interface

  contains

  !> Initializes the module by getting the module parameters from the reader 
  !! object.
  !! @param reader the object which gets parameters by name e.g. from an input 
  !! file.
  !! @param simbox the simulation box defining the boundary conditions and 
  !! volume.
  !! @param particles the particles to be simulated.
  !!
  subroutine initparameterizer(reader, simbox, particles, nbrs)
    type(parameterizer), intent(in) :: reader
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_nbrlist), target, intent(inout), optional :: nbrs
    logical :: overlap
    call energy_init(reader)
    call getparameter(reader, 'scaling_type', scalingtype)
    call parsescalingtype(scalingtype)
    call getparameter(reader, 'temperature', temperature)
    call getparameter(reader, 'pressure', pressure)
    call getparameter(reader, 'move_ratio', moveratio)
    call getparameter(reader, 'scaling_ratio', scalingratio)
    call getparameter(reader, 'largest_translation', largesttranslation)
    call getparameter(reader, 'smallest_translation', smallesttranslation)
    call getparameter(reader, 'largest_rotation', largestrotation)
    call getparameter(reader, 'smallest_rotation', smallestrotation)
    call getparameter(reader, 'max_scaling', maxscaling)
    call getparameter(reader, 'pt_ratio', ptratio)
    call getparameter(reader, 'vol_ratio', volratio)
    call getparameter(reader, 'radial_ratio', radialratio)
    call getparameter(reader, 'nmovetrials', nmovetrials)
    call getparameter(reader, 'nscalingtrials', nscalingtrials)
    call getparameter(reader, 'nacceptedscalings', nacceptedscalings)
    !! The initializations below may affect restart so that it does not result
    !! in the same simulation. 
    nacceptedmoves = 0
    nacceptedscalings = 0
    currentvolume = volume(simbox)
    if (volratio < 1) volratio = size(particles)   
    allocate(nbrlist)
    if (present(nbrs)) then
      nbrlist=>nbrs
    else
      call pnl_init(reader)
      nbrlist = create_nbrlist(simbox, particles)
    end if
    call update(nbrlist, simbox, particles)
    call potentialenergy(simbox, particles, nbrlist, etotal, overlap)
    if(overlap) then 
      stop 'Overlap when initializing mc_sweep! Stopping.'
    end if
    if (ptratio > 0) call pt_init(reader)
  end subroutine 

  !> Checks and parses the scalingtype string to the scalingtypes array. 
  !!
  !! @param scalingtype the string to be parsed.
  !! 
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
  end subroutine

  !> Returns the scalingtypes array.
  !! 
  function getscalingtypes()
    character(len=200), dimension(:), pointer :: getscalingtypes
    getscalingtypes = scalingtypes
  end function

  !> Writes the module parameters and observables to a output file using the 
  !! writer object.
  !!
  !! @param writer the object which defines how parameters are written in a 
  !! file.
  !! 
  subroutine mc_sweep_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'MC sweep parameters')
    call writeparameter(writer, 'move_ratio', moveratio)
    call writeparameter(writer, 'scaling_ratio', scalingratio)
    call writeparameter(writer, 'largest_translation', largesttranslation)
    call writeparameter(writer, 'smallest_translation', smallesttranslation)
    call writeparameter(writer, 'largest_rotation', largestrotation)
    call writeparameter(writer, 'smallest_rotation', smallestrotation)
    call writeparameter(writer, 'max_scaling', maxscaling)
    call writeparameter(writer, 'pt_ratio', ptratio)
    call writeparameter(writer, 'vol_ratio', volratio)
    call writeparameter(writer, 'radial_ratio', radialratio)
    call writeparameter(writer, 'nacceptedradial', nacceptedradial)
    call writeparameter(writer, 'nradialtrials', nradialtrials)
    call writeparameter(writer, 'scaling_type', scalingtype)
    call writeparameter(writer, 'nmovetrials', nmovetrials)
    call writeparameter(writer, 'nacceptedmoves', nacceptedmoves)
    call writeparameter(writer, 'current_move_ratio', &
    real(nacceptedmoves, dp)/real(nmovetrials, dp))
    call writeparameter(writer, 'nscalingtrials', nscalingtrials)
    call writeparameter(writer, 'nacceptedscalings', nacceptedscalings)
    call writeparameter(writer, 'current_scaling_ratio', &
    real(nacceptedscalings, dp)/real(nscalingtrials, dp))
    call writeparameter(writer, 'adjusttype', adjusttype)
    call writeparameter(writer, 'pressure', pressure)
    call writeparameter(writer, 'temperature', temperature)
    call writeparameter(writer, 'volume', currentvolume)
    call writeparameter(writer, 'enthalpy', etotal + currentvolume * pressure)
    call writeparameter(writer, 'total_energy', etotal)
    call pnl_writeparameters(writer)
    call energy_writeparameters(writer)
    call pt_writeparameters(writer)
!    call genvoltrial_writeparameters(writer)
  end subroutine


  !> Runs one sweep of Metropolis Monte Carlo updates to the system. This
  !! Parallel tempering NPT-ensemble sweep consists of trials of particle 
  !! moves, scaling of the simulation box and a exchange of particle system 
  !! coordinates with the particle system in the neighbouring temperature.
  !! 
  !! @param simbox the simulation box defining the limits of particle 
  !! coordinates.
  !! @param particles the array of particles.
  !!  
  subroutine sweep(simbox, particles, genstates, isweep)    
    type(particledat), dimension(:), intent(inout) :: particles
    type(poly_box), intent(inout) :: simbox
    type(rngstate), intent(inout) :: genstates(0:)
    integer, intent(in) :: isweep
    integer :: i
    logical :: overlap
    integer :: irss
    integer :: ivolmove
    !! If the volratio parameter is not given, make volume move once a sweep.
    if (volratio < 1) volratio = size(particles) ! == once per sweep
    if (radialratio < 1) radialratio = 2*size(particles) ! == never
    !if (ptratio < 1) ptratio = size(particles) ! == once per sweep
    !! Particle selected for Random Sequential Skipping:
    irss = int(rng(genstates(0)) * real(size(particles), dp)) + 1
    call update(nbrlist, simbox, particles)
    call potentialenergy(simbox, particles, nbrlist, etotal, overlap)
    if(overlap) then 
      stop 'Overlap when entering sweep! Stopping.'
    end if
    call make_particle_moves(simbox, particles, genstates)
    do ivolmove = 1, size(scalingtypes)
      call movevol(simbox, particles, scalingtypes(ivolmove), genstates(0))
    end do
    !if (mod(i, radialratio) == 0) then
    !  call radialscaling(simbox, particles, genstate)
    !end if
    if (mod(isweep, ptratio) == 0) then
      call makeptmove(simbox, particles, genstates(0))
    end if 
    currentvolume = volume(simbox)
  end subroutine 

  subroutine make_particle_moves(simbox, particles, genstates)
    type(poly_box), intent(in) :: simbox
    type(particledat), intent(inout) :: particles(:)
    type(rngstate), intent(inout) :: genstates(0:)
    integer, allocatable :: indices(:, :)
    integer :: i, thread_id = 0
    !$ integer :: j, ix, iy, iz
    !$ type(simplelist), pointer :: sl

    !! Use domain decomposition if OpenMP parallelization is available:
    !$ if (associated(nbrlist%sl)) then
      !! :TODO: Find out if parallel random number generation is a problem.
      !! :TODO: put OpenMP pragmas here to parallelize the loop below.
      sl=>nbrlist%sl
      !! The allocation below is needed only when the nbrlist has been updated.
      !! Could be optimized.
      !!$ if (allocated(indices)) deallocate(indices)
      !$ allocate(indices(maxval(sl%counts), max(sl%nx/2, 1) * max(sl%ny/2, 1) * max(sl%nz/2, 1))) 

      !! Loop over cells. This can be thought of as looping through a 2 x 2 x 2 cube
      !! of cells.
      !$ do ix=0, min(1, sl%nx-1)
      !$ do iy=0, min(1, sl%ny-1)
      !$ do iz=0, min(1, sl%nz-1)
        !$ indices = reshape(sl%indices(:, ix:sl%nx-1:2, iy:sl%ny-1:2, iz:sl%nz-1:2), (/size(indices(:,1)), size(indices(1,:))/))
        !$OMP PARALLEL
        !$ thread_id = omp_get_thread_num()
        !$ if (size(genstates) /= omp_get_num_threads()) write(*,*) size(genstates), omp_get_num_threads()
        !$OMP DO 
        !$ do i = 1, size(indices(1,:))
          !!write(*, *) "thread:", thread_id, ix, iy, iz, pack(indices(:,i), indices(:,i) > 0)
          !$ do j = 1, size(pack(indices(:,i), indices(:,i) > 0))
            !$ call moveparticle(simbox, particles, indices(j,i), genstates(thread_id))
          !$ end do
        !$ end do
        !$OMP END DO
        !$OMP END PARALLEL
      !$ end do 
      !$ end do
      !$ end do

    !$ else 
      !! Use regular looping if no OpenMP or no cell list.
      do i = 1, size(particles)
        call moveparticle(simbox, particles, i, genstates(thread_id))
      end do
    !$ end if
  end subroutine

  !> Makes one parallel tempering trial move. This is a convenience routine
  !! to enhance readability of code. 
  !!
  !! @param simbox the simulation box.
  !! @param particles the array of particles.
  !!
  subroutine makeptmove(simbox, particles, genstate)
    type(poly_box), intent(inout) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    type(rngstate), intent(inout) :: genstate
    real(dp) :: beta
    real(dp) :: enthalpy
    integer :: nparticles
    logical :: isaccepted
    isaccepted=.false.
    beta = 1._dp/temperature
    enthalpy = etotal + pressure * volume(simbox)
    nparticles = size(particles)
    call pt_move(beta, enthalpy, particles, nparticles, simbox, genstate, isaccepted)
    if (isaccepted) then
      etotal = enthalpy - pressure * volume(simbox)
      !! One way to get rid of explicit list update would be to use some kind 
      !! of observing system between the particlearray and the neighbour list. 
      currentvolume = volume(simbox)
      call update(nbrlist, simbox, particles)
    end if
  end subroutine

  !> Performs a trial move of particles(i).
  !!
  !! @param simbox the simulation box where the particle resides
  !! @param particles the array of particles which particle i is part of and
  !! which contains all the particles that interact with particle i.
  !! @param i the index of particle to be moved.
  !! 
  subroutine moveparticle(simbox, particles, i, genstate)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(in) :: i
    type(rngstate), intent(inout) :: genstate
    type(particledat) :: newparticle
    type(particledat) :: oldparticle
    logical :: overlap
    real(dp) :: enew
    real(dp) :: eold
    enew = 0._dp
    eold = 0._dp
    overlap = .false.
    !! :TODO: Change moving of particles to a separate object which 
    !! :TODO: gets the whole particle array (and possibly the box too).
    newparticle = particles(i)
    call move(newparticle, genstate)
    call setposition(newparticle, minimage(simbox, position(newparticle)))
    oldparticle = particles(i) 
    particles(i) = newparticle
    call potentialenergy(simbox, particles, nbrlist, i, enew, overlap)
    particles(i) = oldparticle
    if(.not. overlap) then 
      call potentialenergy(simbox, particles, nbrlist, i, eold, overlap)
      if (overlap) then
        stop 'moveparticle: overlap with old particle'
      end if
      if(acceptchange(eold, enew, genstate)) then
        particles(i) = newparticle
        etotal = etotal + (enew - eold)
        nacceptedmoves = nacceptedmoves + 1
        call update(nbrlist, simbox, particles, i)
      end if
    end if 
    nmovetrials = nmovetrials + 1
  end subroutine

  !> Returns the simulation temperature.
  !! 
  pure function gettemperature() result(temp)
    real(dp) :: temp
    temp = temperature
  end function

  !> Sets the simulation temperature.
  !!  
  !! @param temperaturein the new temperature
  !!
  subroutine settemperature(temperaturein)
    real(dp), intent(in) :: temperaturein
    temperature = temperaturein
  end subroutine

  !> Performs a trial volume scaling which scales the @param simbox and all the
  !! positions of @param particles. For more information see for example Allen
  !! and Tildesley: Computer Simulation of Liquids the chapter about NPT 
  !! ensemble Monte Carlo.
  !! 
  !! @param simbox the simulation box.
  !! @param particles the particles in the simulation box.
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
    nparticles = size(particles)
    overlap = .false.
    Vo = volume(simbox)
    oldbox = simbox
    scaling = genvoltrial_scale(simbox, maxscaling, genstate, trim(adjustl(scalingtype)))
    call scalepositions(oldbox, simbox, particles, nparticles) 
    Vn = volume(simbox)
    call potentialenergy(simbox, particles, nbrlist, totenew, overlap)    
    call scalepositions(simbox, oldbox, particles, nparticles)
    if (overlap) then
      simbox = oldbox  
    else
      boltzmannn = totenew + pressure * Vn - real(nparticles, dp) * &
        temperature * log(Vn)  
      boltzmanno = etotal + pressure * Vo - real(nparticles, dp) * &
        temperature * log(Vo)
      isaccepted = acceptchange(boltzmanno, boltzmannn, genstate)
      if (isaccepted) then
        !! Scale back to new configuration
        call scalepositions(oldbox, simbox, particles, nparticles)
        etotal = totenew
        nacceptedscalings = nacceptedscalings + 1
        currentvolume = volume(simbox)
        call update(nbrlist, simbox, particles)
      else 
        simbox = oldbox
      end if 
    end if
    nscalingtrials = nscalingtrials + 1
  end subroutine

  !! Makes a trial scaling of @p particles' radial coordinates. Does not scale
  !! @p simbox dimensions. 
  !! 
  !! @p simbox the box that contains the particles
  !! @p particles the array of particles in simbox
  !! @p genstate the random number generator state
  !! 
  subroutine radialscaling(simbox, particles, genstate)    
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    type(rngstate), intent(inout) :: genstate
    integer :: nparticles 
    logical :: overlap
    logical :: isaccepted
    real(dp) :: totenew
    type(poly_box) :: tempbox
    real(dp), dimension(3) :: scaling
    nparticles = size(particles)
    overlap = .false.
    tempbox = simbox
    scaling = genvoltrial_scale(tempbox, maxscaling, genstate, 'xy')
    call scalepositions(simbox, tempbox, particles, nparticles)
    call potentialenergy(simbox, particles, nbrlist, totenew, overlap)    
    !! Scale back to old coordinates
    call scalepositions(tempbox, simbox, particles, nparticles)
    if (.not. overlap) then
      isaccepted = acceptchange(etotal, totenew, genstate)
      if (isaccepted) then
        !! Scale back to new configuration
        call scalepositions(simbox, tempbox, particles, nparticles)
        etotal = totenew
        nacceptedradial = nacceptedradial + 1
        call update(nbrlist, simbox, particles)
      end if 
    end if
    nradialtrials = nradialtrials + 1
  end subroutine

  !> Returns a boolean value that tells if the Metropolis trial move with the
  !! given old and new energies is accepted. .true. means yes.
  !!
  !! @param oldenergy the energy of the system before the move.
  !! @param newenergy the energy of the system after the trial move.
  !!
  function acceptchange(oldenergy, newenergy, genstate) &
    result(isaccepted)    
    logical :: isaccepted
    real(dp), intent(in) :: oldenergy
    real(dp), intent(in) :: newenergy
    type(rngstate), intent(inout) :: genstate
    real(dp) :: dE
    dE = newenergy - oldenergy
    if(dE < 0._dp) then
      isaccepted = .true.
    else
      isaccepted = (rng(genstate) < exp(-dE/temperature))
    end if  
  end function

  !> Adjusts the maximum values for trial moves of particles and trial scalings
  !! of the simulation volume. Should be used only during equilibration run.
  !! 
  subroutine updatemaxvalues
    real(dp) :: newdthetamax
    real(dp) :: newdximax
    !! Adjust scaling
    maxscaling = newmaxvalue(nscalingtrials, nacceptedscalings, scalingratio,&
    maxscaling)
    call getmaxmoves(newdximax, newdthetamax)
    !! Adjust translation
    newdximax = newmaxvalue(nmovetrials, nacceptedmoves, moveratio, newdximax)
    if (newdximax > largesttranslation) then
      !write(*, *) 'Tried to increase maxtrans above largesttranslation.'
      newdximax = largesttranslation
    else if(newdximax < smallesttranslation) then
      !write(*, *) 'Tried to decrease maxtrans below smallesttranslation.'
      newdximax = smallesttranslation
    end if
    !! Adjust rotation
    newdthetamax = newmaxvalue(nmovetrials, nacceptedmoves, moveratio, &
    newdthetamax)
    if (newdthetamax > largestrotation) then 
      newdthetamax = largestrotation
    else if (newdthetamax < smallestrotation) then
      newdthetamax = smallestrotation
    end if
    call setmaxmoves(newdximax, newdthetamax)
  end subroutine
  
  !> Returns a new trial move parameter value calculated from the desired 
  !! acceptance ratio
  !! 
  !! @param ntrials the total number of trials of the kind of move in question.
  !! @param naccepted number of accepted trials.
  !! @param desiredratio the desired acceptance ratio for trial moves.
  !! @param oldvalue the old value of the parameter setting the maximum size 
  !! for the trial move in question.
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
  end function

  !> Returns the simulation pressure in reduced units.
  !!
  !! @return the simulation pressure. 
  !! 
  function getpressure()
    real(dp) :: getpressure
    getpressure = pressure
  end function

  !> Resets the counters that are used to monitor acceptances of trial moves. 
  !!
  subroutine resetcounters
    nacceptedmoves = 0
    nmovetrials = 0
    nscalingtrials = 0 
    nacceptedscalings = 0
    call pt_resetcounters
  end subroutine

  !> Scales the positions of @param particles with the same factors that are 
  !! used to scale the simulation box dimensions from oldbox to newbox. To be 
  !! used with NPT ensemble Metropolis Monte Carlo.
  !! 
  !! @param oldbox the simulation box before scaling.
  !! @param newbox the simulation box after scaling.
  !! @param particles the particles for which the positions are scaled.
  !! @param nparticles the number of particles. 
  !!
  subroutine scalepositions(oldbox, newbox, particles, nparticles)
    implicit none
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
  end subroutine

end module
