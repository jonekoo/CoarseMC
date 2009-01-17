module mc_sweep
  use nrtype, only : dp
  use energy, only: total_energy, potential_energy
  use cylinder, only: setLz, setR, scale, volume, getHeight
  use particle, only : particledat, move, setmaxmoves, getmaxmoves
  use mtmod, only : grnd
  use verlet, only : updatelist, totpairV



  public :: init
  public :: sweep
  public :: save_state
  public :: load_state
  public :: updatemaxvalues



  private 

  integer, save :: n_accepted_moves_
  integer, save :: n_accepted_scalings_
  real(dp), save :: max_scaling_
  integer, save :: volume_change_type_ 
  real(dp), save :: temperature_
  real(dp), save :: pressure_
  namelist /mcsweep_nml/ temperature_, pressure_, n_accepted_moves_, &
    & n_accepted_scalings_, max_scaling_, volume_change_type_
 


  contains

  subroutine init(volume_change_type, temperature, pressure)
    implicit none
    integer, intent(in) :: volume_change_type
    real(dp), intent(in) :: temperature
    real(dp), intent(in) :: pressure
    volume_change_type_ = volume_change_type
    temperature_ = temperature
    pressure_ = pressure 
    max_scaling_ = 0.1
    n_accepted_moves_ = 0
    n_accepted_scalings_ = 0
  end subroutine init



  subroutine save_state(write_unit)
    implicit none 
    integer, intent(in) :: write_unit
    write(write_unit, NML = mcsweep_nml)
  end subroutine save_state



  subroutine load_state(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    read(read_unit, NML = mcsweep_nml)
  end subroutine



  subroutine sweep(particles, n_particles)
    implicit none
    intrinsic log, real
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(in) :: n_particles
    integer :: i
    real(dp) :: Eold = 0.0, Enew = 0.0, totEold = 0.0, totEnew = 0.0
    real(dp) :: Lz0, Lzn
    real(dp) :: V0, Vn
    logical :: accept = .true.
    logical :: overlap
    type(particledat) :: newparticle 
    !! Trial moves of particles
    do i = 1, n_particles
      Enew = 0
      Eold = 0
      overlap = .false.
      call move(particles(i), newparticle)
      call potential_energy(particles, n_particles, newparticle, i, Enew, & 
        overlap)
      if(.not. overlap) then 
        call potential_energy(particles, n_particles, particles(i), i, & 
          Eold, overlap)   
        accept = acceptchange(Eold, Enew, overlap)       
        if(accept) then
          particles(i) = newparticle
          n_accepted_moves_ = n_accepted_moves_ + 1
        end if
      end if 
    end do
    !! Trial volume change
    overlap = .false.
    V0 = volume()
    Lz0 = getHeight()
    Lzn = Lz0 + (2.0 * grnd() - 1.0) * max_scaling_
    call changeLz(particles, n_particles, Lz0, Lzn);
    Vn = volume()
    call totpairV(particles, n_particles, totEnew, overlap)
    call changeLz(particles, n_particles, Lzn, Lz0);
    if (overlap) return; !! If the new state results in an overlap, 
                         !! old state is restored
    call totpairV(particles, n_particles, totEold, overlap)  
    totEnew = totEnew+pressure_*Vn-real(n_particles)*temperature_*log(Vn)   
    totEold = totEold+pressure_*V0-real(n_particles)*temperature_*log(V0)
    accept = acceptchange(totEold, totEnew, overlap)
    if (accept) then
      call changeLz(particles, n_particles, Lz0, Lzn);
      n_accepted_scalings_ = n_accepted_scalings_ + 1
    end if 
    call updatelist(particles, n_particles)
  end subroutine sweep



  !! Funktio, joka uuden ja vanhan energian perusteella
  !! päättää, hyväksytäänkö muutos
  logical function acceptchange(oldenergy, newenergy, overlap) &
    result(is_accepted)
    implicit none
    intrinsic exp
    real(dp), intent(in) :: oldenergy
    real(dp), intent(in) :: newenergy
    real(dp) :: u, dE
    logical, intent(in) :: overlap    
    if(overlap) then 
      is_accepted=.false.
      return;
    end if
    dE = newenergy - oldenergy
    is_accepted = .false.
    if(dE < 0) then
      is_accepted = .true.
    else
      u = grnd()
      if ( u < exp(-dE/temperature_) ) then 
        is_accepted = .true.
      end if
    end if  
  end function acceptchange



  !! Palauttaa hyväksymisuhteet siirron ja kierron yhdistelmälle,
  !! sekä tilavuuden muutokselle
  subroutine ratios(Nparticles,period,movratio,volratio)
    implicit none    
    real, intent(out) :: movratio, volratio
    integer, intent(in) :: Nparticles, period
    movratio = real(n_accepted_moves_)/real(Nparticles*period)
    volratio = real(n_accepted_scalings_)/real(period)   
  end subroutine ratios



  !! Päivittää maksimimuutosarvot koordinaateille/kulmille ja 
  !! sylinterin säteelle.
  subroutine updatemaxvalues(Nparticles, period)
    implicit none
    real :: mratio, vratio
    logical :: volinc, movinc
    real(dp) :: olddximax, olddthetamax
    real(dp) :: newdximax, newdthetamax
    integer, intent(in) :: Nparticles, period
    call ratios(Nparticles, period, mratio, vratio)      
    !Jos tilavuuden muutoksista on hyväksytty yli 25%,
    !kasvatetaan säteen maksimimuutosarvoa. Vastaavasti
    !siirrolle/kierrolle rajana 33%.
    volinc = (0.25 < vratio)
    movinc = (0.33 < mratio)
    !Nollataan hyväksyntälaskurit
    n_accepted_scalings_ = 0
    n_accepted_moves_ = 0
    max_scaling_ = newmaxvalue(volinc, max_scaling_)      
    call getmaxmoves(olddximax, olddthetamax)
    newdthetamax = newmaxvalue(movinc, olddthetamax)
    newdximax = newmaxvalue(movinc, olddximax)
    call setmaxmoves(newdximax, newdthetamax)
  end subroutine updatemaxvalues
  


  function newmaxvalue(increase, oldvalue) result(newvalue)
    implicit none
    logical, intent(in) :: increase
    real(dp), intent(in) :: oldvalue
    real(dp) :: newvalue
    real(dp), parameter :: multiplier = 1.05
    if (increase) then
      newvalue = oldvalue * multiplier
    else 
      newvalue = oldvalue / multiplier
    end if
  end function newmaxvalue
    


  subroutine changeR(particlearray, oldR,newR)
    implicit none
    real(dp),intent(in) :: oldR,newR
    type(particledat), dimension(:), pointer :: particlearray 
    integer :: i
    real(dp) :: x0,y0        
    call setR(newR)
    do i = 1, size(particlearray)
      x0 = particlearray(i)%x
      y0 = particlearray(i)%y
      call scale(x0, y0, particlearray(i)%x, particlearray(i)%y, oldR, newR)
    end do
  end subroutine changeR



  subroutine changeLz(particles, n_particles, oldLz, newLz)
    implicit none
    type(particledat),dimension(:), intent(inout) :: particles
    integer, intent(in) :: n_particles
    real(dp), intent(in) :: oldLz, newLz
    integer :: i
    call setLz(newLz)
    do i = 1, n_particles
      particles(i)%z=newLz/oldLz*particles(i)%z
    end do
  end subroutine changeLz



end module mc_sweep
