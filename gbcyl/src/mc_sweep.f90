module mc_sweep
  use nrtype
  use energy
  use class_poly_box
  use particle
  use mtmod
!  use verlet
  use cell
  use cell_energy, only: new_list, pairinteractions  
  use pt
  use class_parameterizer
  use class_parameter_writer
  implicit none  
  private 

  public :: init
  public :: sweep
  public :: updatemaxvalues
  public :: getpressure
  public :: mc_sweep_writeparameters

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
    !! The initializations below may affect restart so that it does not result
    !! in the same simulation. 
    nacceptedmoves = 0
    nacceptedscalings = 0
  end subroutine

  subroutine mc_sweep_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'MC sweep parameters')
    call writeparameter(writer, 'volume_change_type', volumechangetype)
    call writeparameter(writer, 'pt_low', ptlow)
    call writeparameter(writer, 'pt_high', pthigh)
    call writeparameter(writer, 'pressure', pressure)
    call writeparameter(writer, 'move_ratio', moveratio)
    call writeparameter(writer, 'scaling_ratio', scalingratio)
    call writeparameter(writer, 'largest_translation', largesttranslation)
    call writeparameter(writer, 'max_scaling', maxscaling)
    call writeparameter(writer, 'pt_ratio', ptratio)
  end subroutine

  subroutine sweep(simbox, particles, nbrlist)
    type(particledat), dimension(:), intent(inout) :: particles
    type(poly_box), intent(inout) :: simbox
    type(list), intent(inout) :: nbrlist
    integer :: i
    logical :: overlap
    integer :: ivolmove
    real(dp) :: etotalold
    !! Trial moves of particles and one particle move replaced with volume move
    ivolmove = int(grnd() * real(size(particles), dp)) + 1
    etotalold = etotal   
    nbrlist = new_list(simbox, particles)
    call totalenergy(simbox, particles, nbrlist, etotal, overlap)
    if(overlap) then 
      stop 'Overlap when entering sweep! Stopping.'
    end if
    do i = 1, size(particles)
      if (i == ivolmove) then
        !! Replace random particle update with a volume update.
        call movevol(simbox, particles, nbrlist)
      else
        call moveparticle(simbox, particles, nbrlist, i)
      end if
      if (mod(i, ptratio) == 0) then
        call makeptmove(simbox, particles, nbrlist)
      end if 
    end do
  end subroutine

  subroutine makeptmove(simbox, particles, nbrlist)
    type(poly_box), intent(inout) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    type(list), intent(inout) :: nbrlist
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
    nbrlist = new_list(simbox, particles(1:size(particles)))
  end subroutine

  subroutine moveparticle(simbox, particles, nbrlist, i)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    type(list), intent(inout) :: nbrlist
    integer, intent(in) :: i
    type(particledat) :: newparticle
    type(particledat) :: oldparticle
    logical :: overlap
    real(dp) :: enew
    real(dp) :: eold
    logical :: isaccepted
    enew = 0._dp
    eold = 0._dp
    overlap = .false.
    !! :TODO: Change moving of particles for a separate object which 
    !! :TODO: gets the whole particle array.
    call move(particles(i), newparticle)
    call setposition(newparticle, minimage(simbox, &
    (/0._dp, 0._dp, 0._dp/), position(newparticle)))
    oldparticle = particles(i) 
    particles(i) = newparticle
    call potentialenergy(simbox, particles, nbrlist, i, enew, overlap)
    if(.not. overlap) then 
      particles(i) = oldparticle
      call potentialenergy(simbox, particles, nbrlist, i, eold, overlap)
      if (overlap) then
        stop 'sweep: overlap with old particle'
      end if
      isaccepted = acceptchange(eold, enew)       
      if(isaccepted) then
        particles(i) = newparticle
        etotal = etotal + (enew - eold)
        nacceptedmoves = nacceptedmoves + 1
      end if
    end if 
  end subroutine

  !subroutine moveparticleinsystem(simbox, particles, i)
  !  type(poly_box), intent(in) :: simbox
  !  type(particledat), dimension(:), intent(inout) :: particles
  !  integer, intent(in) :: i
  !  type(particledat) :: newparticle
  !  type(particledat) :: oldparticle
  !  logical :: overlap
  !  real(dp) :: enew
  !  real(dp) :: eold
  !  logical :: isaccepted
  !  enew = 0._dp
  !  eold = 0._dp
  !  overlap = .false.
    !! :TODO: Change moving of particles for a separate object which 
    !! :TODO: gets the whole particle array.
    !call moveparticle(particleiterator) 
    !call potentialenergy(particleiterator, enew, overlap) 
    !if(.not. overlap) then 
    !  call undomove(particleiterator)
    !  call potentialenergy(particleiterator, eold, overlap)
    !  if (overlap) then
    !    stop 'sweep: overlap with old particle'
    !  end if
    !  isaccepted = acceptchange(eold, enew)       
    !  if (isaccepted) then
    !    call redomove(particleiterator)
    !    etotal = etotal + (enew - eold) !! :TODO: Remove this.
    !    nacceptedmoves = nacceptedmoves + 1
    !  end if
    !end if 
  !end subroutine

  subroutine movevol(simbox, particles, nbrlist)
    type(poly_box), intent(inout) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    type(list), intent(in) :: nbrlist
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
    call totalenergy(simbox, particles, nbrlist, totenew, overlap)
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

  !! Palauttaa hyv‰ksymisuhteet siirron ja kierron yhdistelm‰lle,
  !! sek‰ tilavuuden muutokselle
  !! 
  subroutine ratios(Nparticles, period, movratio, volratio)
    real(dp), intent(out) :: movratio, volratio
    integer, intent(in) :: Nparticles, period
    movratio = real(nacceptedmoves, dp)/real(Nparticles*period, dp)
    volratio = real(nacceptedscalings, dp)/real(period, dp)   
  end subroutine

  !! P‰ivitt‰‰ maksimimuutosarvot koordinaateille/kulmille ja 
  !! sylinterin s‰teelle.
  !!
  subroutine updatemaxvalues(Nparticles, period)
    real(dp) :: mratio, vratio
    real(dp) :: olddximax, olddthetamax
    real(dp) :: newdximax, newdthetamax
    integer, intent(in) :: Nparticles, period
    !real(dp), save :: lastmratio = 0._dp
    !real(dp), save :: dmacceptance = 0._dp
    !real(dp) :: dtransrot
    call ratios(Nparticles, period, mratio, vratio)      
    !! Jos tilavuuden muutoksista on hyv‰ksytty yli 25%,
    !! kasvatetaan s‰teen maksimimuutosarvoa. Vastaavasti
    !! siirrolle/kierrolle rajana 33%.
    !! Nollataan hyv‰ksynt‰laskurit
    nacceptedscalings = 0
    nacceptedmoves = 0
    maxscaling = newmaxvalue(vratio > scalingratio, maxscaling)
    call getmaxmoves(olddximax, olddthetamax)
    !if(adjusttype == 1 .or. adjusttype == 2) then
      newdthetamax = newmaxvalue(mratio > moveratio, olddthetamax)
      newdximax = newmaxvalue(mratio > moveratio, olddximax)
    !end if
    !if (adjusttype == 2) then
      !! adjust ratio maxtranslation/maxrotation
    !  dmacceptance = mratio - lastmratio
    !  if (dmacceptance*dtransrot > 0._dp) then
    !    newdximax = 1.05_dp*newdximax
    !    dtransrot = 1._dp
    !  else
    !    newdximax = newdximax/1.05_dp
    !    dtransrot = -1._dp
    !  end if
    !  lastmratio = mratio
    !end if
    !! :TODO: write to log if trying to increase max move too much
    !! :TODO: add also smallest possible newdximax? 
    !! :TODO: write to log if trying to decrease below smallest newdximax
    if (newdximax > largesttranslation) then
      !write(*, *) 'Tried to increase maxtrans above largesttranslation.'
      newdximax = largesttranslation
      newdthetamax = largestrotation
    else if(newdximax < smallesttranslation) then
      !write(*, *) 'Tried to decrease maxtrans below smallesttranslation.'
      newdximax = smallesttranslation
      newdthetamax = smallestrotation
    end if
    call setmaxmoves(newdximax, newdthetamax)
  end subroutine
  
  function newmaxvalue(increase, oldvalue) result(newvalue)
    logical, intent(in) :: increase
    real(dp), intent(in) :: oldvalue
    real(dp) :: newvalue
    real(dp), parameter :: multiplier = 1.05_dp
    if (increase) then
      newvalue = oldvalue * multiplier
    else 
      newvalue = oldvalue / multiplier
    end if
  end function

  function getpressure()
    real(dp) :: getpressure
    getpressure = pressure
  end function

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
