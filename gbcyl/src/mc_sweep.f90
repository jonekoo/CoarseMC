module mc_sweep
  use nrtype
  use energy
  use class_poly_box
  use particle
  use mtmod
  use class_poly_nbrlist
  use pt
  use class_parameterizer
  use class_parameter_writer
  use mpi
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

  integer, save :: nacceptedmoves
  integer, save :: nacceptedscalings
  real(dp), save :: maxscaling
  integer, save :: volumechangetype 
  real(dp), save :: temperature
  real(dp), save :: pressure
  integer, save :: adjusttype
  real(dp), save :: moveratio
  real(dp), save :: scalingratio
  real(dp), save :: largesttranslation = 1._dp    !! molecule diameter
  real(dp), save :: smallesttranslation = 0.1_dp !! 1/10 molecule diameter
  real(dp), save :: largestrotation = 1.57_dp     !! ~ pi/4
  real(dp), save :: smallestrotation = 0.157_dp  !! 1/10 largest rotation
  real(dp), save :: ptlow
  real(dp), save :: pthigh
  integer, save :: ptratio = 10
  real(dp), save :: etotal = 0._dp
  real(dp), save :: currentvolume = 0._dp
  integer, save :: nmovetrials = 0
  integer, save :: nscalingtrials = 0

  interface init
    module procedure initparameterizer
  end interface

  contains

  subroutine initparameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'volume_change_type', volumechangetype)
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
    call getparameter(reader, 'nmovetrials', nmovetrials)
    call getparameter(reader, 'nscalingtrials', nscalingtrials)
    call getparameter(reader, 'nacceptedscalings', nacceptedscalings)
    !! The initializations below may affect restart so that it does not result
    !! in the same simulation. 
    nacceptedmoves = 0
    nacceptedscalings = 0
  end subroutine

  subroutine mc_sweep_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'MC sweep parameters')
    call writeparameter(writer, 'volume_change_type', volumechangetype)
    call writeparameter(writer, 'move_ratio', moveratio)
    call writeparameter(writer, 'scaling_ratio', scalingratio)
    call writeparameter(writer, 'largest_translation', largesttranslation)
    call writeparameter(writer, 'max_scaling', maxscaling)
    call writeparameter(writer, 'pt_ratio', ptratio)
    call writeparameter(writer, 'nmovetrials', nmovetrials)
    call writeparameter(writer, 'nacceptedmoves', nacceptedmoves)
    call writeparameter(writer, 'current_move_ratio', real(nacceptedmoves, dp)/real(nmovetrials))
    call writeparameter(writer, 'nscalingtrials', nscalingtrials)
    call writeparameter(writer, 'nacceptedscalings', nacceptedscalings)
    call writeparameter(writer, 'current_scaling_ratio', real(nacceptedscalings, dp)/real(nscalingtrials))
    call writeparameter(writer, 'adjusttype', adjusttype)
    call writeparameter(writer, 'smallesttranslation', smallesttranslation)
    call writeparameter(writer, 'largestrotation', largestrotation)
    call writeparameter(writer, 'smallestrotation', smallestrotation)
    call writeparameter(writer, 'pressure', pressure)
    call writeparameter(writer, 'temperature', temperature)
    call writeparameter(writer, 'volume', currentvolume)
    call writeparameter(writer, 'enthalpy', etotal + currentvolume * pressure)
    call writeparameter(writer, 'etotal', etotal)
    !call writeparameter(writer, 'measured_move_ratio', mratio) 
    !call writeparameter(writer, 'measured_scaling_ratio', vratio) 
  end subroutine

  subroutine sweep(simbox, particles, nbrlist)
    type(particledat), dimension(:), intent(inout) :: particles
    type(poly_box), intent(inout) :: simbox
    type(poly_nbrlist), intent(inout) :: nbrlist
    integer :: i
    logical :: overlap
    integer :: ivolmove
    !! Trial moves of particles and one particle move replaced with volume move
    ivolmove = int(grnd() * real(size(particles), dp)) + 1
    call update(nbrlist, simbox, particles)
    call potentialenergy(simbox, particles, nbrlist, etotal, overlap)
    if(overlap) then 
      stop 'Overlap when entering sweep! Stopping.'
    end if
    do i = 1, size(particles)
      if (i == ivolmove) then
        call movevol(simbox, particles, nbrlist)
        call update(nbrlist, simbox, particles)
      else
        call moveparticle(simbox, particles, nbrlist, i)
      end if
      if (mod(i, ptratio) == 0) then
        call makeptmove(simbox, particles, nbrlist)
        call update(nbrlist, simbox, particles)
      end if 
    end do
    currentvolume = volume(simbox)
  end subroutine

  subroutine makeptmove(simbox, particles, nbrlist)
    type(poly_box), intent(inout) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    type(poly_nbrlist), intent(inout) :: nbrlist
    real(dp) :: beta
    real(dp) :: enthalpy
    integer :: nparticles
    beta = 1._dp/temperature
    enthalpy = etotal + pressure * volume(simbox)
    nparticles = size(particles)
    call pt_move(beta, enthalpy, particles, nparticles, simbox)
    etotal = enthalpy - pressure * volume(simbox)
    !! One way to get rid of explicit list update would be to use some kind 
    !! of observing system between the particlearray and the neighbour list. 
  end subroutine

  subroutine moveparticle(simbox, particles, nbrlist, i)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    type(poly_nbrlist), intent(in) :: nbrlist
    integer, intent(in) :: i
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
    call move(newparticle)
    call setposition(newparticle, minimage(simbox, position(newparticle)))
    oldparticle = particles(i) 
    particles(i) = newparticle
    call potentialenergy(simbox, particles, nbrlist, i, enew, overlap)
    particles(i) = oldparticle
    if(.not. overlap) then 
      call potentialenergy(simbox, particles, nbrlist, i, eold, overlap)
      if (overlap) then
        stop 'sweep: overlap with old particle'
      end if
      if(acceptchange(eold, enew)) then
        particles(i) = newparticle
        etotal = etotal + (enew - eold)
        nacceptedmoves = nacceptedmoves + 1
      end if
    end if 
    nmovetrials = nmovetrials + 1
  end subroutine

  pure function gettemperature() result(temp)
    real(dp) :: temp
    temp = temperature
  end function

  subroutine settemperature(temperaturein)
    real(dp), intent(in) :: temperaturein
    temperature = temperaturein
  end subroutine

  subroutine movevol(simbox, particles, nbrlist)
    type(poly_box), intent(inout) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    type(poly_nbrlist), intent(in) :: nbrlist
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
    scaling = scale(simbox, maxscaling, grnd)
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
      isaccepted = acceptchange(boltzmanno, boltzmannn)
      if (isaccepted) then
        !! Scale back to new configuration
        call scalepositions(oldbox, simbox, particles, nparticles)
        etotal = totenew
        nacceptedscalings = nacceptedscalings + 1
      else 
        simbox = oldbox
      end if 
    end if
    nscalingtrials = nscalingtrials + 1
  end subroutine

  !! Funktio, joka uuden ja vanhan energian perusteella
  !! p‰‰tt‰‰, hyv‰ksyt‰‰nkˆ muutos
  !!
  function acceptchange(oldenergy, newenergy) &
    result(isaccepted)
    intrinsic exp
    logical :: isaccepted
    real(dp), intent(in) :: oldenergy
    real(dp), intent(in) :: newenergy
    real(dp) :: dE
    dE = newenergy - oldenergy
    if(dE < 0._dp) then
      isaccepted = .true.
    else
      isaccepted = (grnd() < exp(-dE/temperature))
    end if  
  end function

  !! Adjusts the maximum values for trial moves of particles and trial scalings
  !! of the simulation volume.
  !! 
  subroutine updatemaxvalues
    real(dp) :: newdthetamax
    real(dp) :: newdximax
    !! Adjust scaling
    maxscaling = newmaxvalue(nscalingtrials, nacceptedscalings, scalingratio, maxscaling)

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
    newdthetamax = newmaxvalue(nmovetrials, nacceptedmoves, moveratio, newdthetamax)
    if (newdthetamax > largestrotation) then 
      newdthetamax = largestrotation
    else if (newdthetamax < smallestrotation) then
      newdthetamax = smallestrotation
    end if
    call setmaxmoves(newdximax, newdthetamax)
  end subroutine
  
  function newmaxvalue(ntrials, naccepted, desiredratio, oldvalue) result(newvalue)
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

  function getpressure()
    real(dp) :: getpressure
    getpressure = pressure
  end function

  subroutine resetcounters
    nacceptedmoves = 0
    nmovetrials = 0
    nscalingtrials = 0 
    nacceptedscalings = 0
    call pt_resetcounters
  end subroutine

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

end module mc_sweep
