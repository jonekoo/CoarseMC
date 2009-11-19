module mc_sweep
  use nrtype, only : dp
  use energy, only: total_energy, potential_energy
  use class_poly_box
  use particle, only : particledat, move, setmaxmoves, getmaxmoves, position, &
    set_position
  use mtmod, only : grnd
  use verlet, only : pair_interactions, updatelist
!  use all_pairs, only: pair_interactions
  use pt
  use class_parameterizer
  use class_parameter_writer
  implicit none  
  private 

  public :: init
  public :: sweep
  public :: updatemaxvalues
  public :: pressure
  public :: mc_sweep_write_parameters

  integer, save :: n_accepted_moves_
  integer, save :: n_accepted_scalings_
  real(dp), save :: max_scaling_
  integer, save :: volume_change_type_ 
  real(dp), save :: temperature_
  real(dp), save :: pressure_
  integer, save :: adjust_type_
  real(dp), save :: move_ratio_
  real(dp), save :: scaling_ratio_
  real(dp), save :: largest_translation_ = 1._dp    !! molecule diameter
  real(dp), save :: smallest_translation_ = 0.01_dp !! 1/100 molecule diameter
  real(dp), save :: largest_rotation_ = 1.57_dp     !! ~ pi/4
  real(dp), save :: smallest_rotation_ = 0.0157_dp  !! 1/100 largest rotation
  real(dp), save :: pt_low
  real(dp), save :: pt_high
  integer, save :: pt_ratio_ = 10
  real(dp), save :: e_total_ = 0._dp
 
  interface init
    module procedure init_parameterizer
  end interface

  contains

  subroutine init_parameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call get_parameter(reader, 'volume_change_type', volume_change_type_)
    call get_parameter(reader, 'pt_low', pt_low)
    call get_parameter(reader, 'pt_high', pt_high)
    temperature_ = pt_temperature(pt_low, pt_high)
    call get_parameter(reader, 'pressure', pressure_)
    call get_parameter(reader, 'move_ratio', move_ratio_)
    call get_parameter(reader, 'scaling_ratio', scaling_ratio_)
    call get_parameter(reader, 'largest_translation', largest_translation_)
    call get_parameter(reader, 'smallest_translation', smallest_translation_)
    call get_parameter(reader, 'largest_rotation', largest_rotation_)
    call get_parameter(reader, 'smallest_rotation', smallest_rotation_)
    call get_parameter(reader, 'max_scaling', max_scaling_)
    call get_parameter(reader, 'pt_ratio', pt_ratio_)
    !! The initializations below may affect restart so that it does not result
    !! in the same simulation. 
    n_accepted_moves_ = 0
    n_accepted_scalings_ = 0
  end subroutine

  subroutine mc_sweep_write_parameters(writer)
    type(parameter_writer), intent(in) :: writer
    call write_comment(writer, 'MC sweep parameters')
    call write_parameter(writer, 'volume_change_type', volume_change_type_)
    call write_parameter(writer, 'pt_low', pt_low)
    call write_parameter(writer, 'pt_high', pt_high)
    call write_parameter(writer, 'pressure', pressure_)
    call write_parameter(writer, 'move_ratio', move_ratio_)
    call write_parameter(writer, 'scaling_ratio', scaling_ratio_)
    call write_parameter(writer, 'largest_translation', largest_translation_)
    call write_parameter(writer, 'max_scaling', max_scaling_)
    call write_parameter(writer, 'pt_ratio', pt_ratio_)
  end subroutine

  subroutine sweep(particles, n_particles, simbox)
    intrinsic log, real
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(inout) :: n_particles
    type(poly_box), intent(inout) :: simbox
    integer :: i
    real(dp) :: Eold
    real(dp) :: Enew
    logical :: accept
    logical :: overlap
    type(particledat) :: newparticle 
    real(dp) :: beta
    real(dp) :: enthalpy
    integer :: i_vol_move
    real(dp) :: e_total_old
    !! Trial moves of particles and one particle move replaced with volume move
    i_vol_move = int(grnd() * real(n_particles, dp)) + 1
    e_total_old = e_total_
    call total_energy(particles, n_particles, simbox, e_total_, overlap)
    if(overlap) then 
      stop 'Overlap when entering sweep! Stopping.'
    end if
    !write(*, *) 'These should be the same: ', e_total_old, e_total_
    do i = 1, n_particles
      if (i == i_vol_move) then
        call move_vol(simbox, particles, n_particles)
      else
        Enew = 0._dp
        Eold = 0._dp
        overlap = .false.
        !! :TODO: Change moving of particles for a separate object which gets the whole particle array!
        call move(particles(i), newparticle)
        call set_position(newparticle, min_image(simbox, (/0._dp, 0._dp, 0._dp/), position(newparticle))) 
        call potential_energy(particles, n_particles, newparticle, i, simbox, &
          Enew, overlap)
        if(.not. overlap) then 
          call potential_energy(particles, n_particles, particles(i), i, simbox, &
            Eold, overlap)   
          if (overlap) then
            stop 'sweep: overlap with old particle'
          end if
          accept = acceptchange(Eold, Enew)       
          if(accept) then
            particles(i) = newparticle
            e_total_ = e_total_ + (Enew - Eold)
            n_accepted_moves_ = n_accepted_moves_ + 1
          end if
        end if 
      end if
      if (mod(i, pt_ratio_) == 0) then
        beta = 1._dp/temperature_
        !call total_energy(particles, n_particles, simbox, enthalpy, overlap)
        !write(*, *) "These should be the same: ", e_total_, enthalpy
        !if (overlap) then 
        !  stop 'Overlap before parallel tempering move!'
        !end if
        enthalpy = e_total_ + pressure_ * volume(simbox)
        call pt_move(beta, enthalpy, particles, n_particles, simbox)
        e_total_ = enthalpy - pressure_ * volume(simbox)
        !! One way to get rid of explicit list update would be to use some kind 
        !! of observing system between the particlearray and the neighbour list. 
        call updatelist(particles, n_particles, simbox) 
      end if 
    end do
  end subroutine

  subroutine move_vol(simbox, particles, n_particles)
    type(poly_box), intent(inout) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(in) :: n_particles
    logical :: overlap
    real(dp) :: V_o, V_n
    logical :: accept
!    real(dp) :: totEold
    real(dp) :: totEnew
    real(dp) :: boltzmann_n
    real(dp) :: boltzmann_o
    type(poly_box) :: oldbox
    overlap = .false.
    V_o = volume(simbox)
    oldbox = simbox
    call scale(simbox, max_scaling_, grnd)
    call scale_positions(oldbox, simbox, particles, n_particles) 
    V_n = volume(simbox)
    call total_energy(particles, n_particles, simbox, totEnew, overlap)
    call scale_positions(simbox, oldbox, particles, n_particles)
    if (overlap) then
      simbox = oldbox  
    else
      !call total_energy(particles, n_particles, oldbox, totEold, overlap)
      !if (overlap) then
      !  stop 'move_vol: overlap in old configuration'
      !end if
      !write(*, *) 'These should be the same: ', e_total_, totEold
      !e_total_ = totEold
      boltzmann_n = totEnew + pressure_ * V_n - real(n_particles, dp) * &
        temperature_ * log(V_n)  
      boltzmann_o = e_total_ + pressure_ * V_o - real(n_particles, dp) * &
        temperature_ * log(V_o)
      accept = acceptchange(boltzmann_o, boltzmann_n)
      if (accept) then
        !! Scale back to new configuration
        !write(*, *) 'Accepted volume move' 
        call scale_positions(oldbox, simbox, particles, n_particles)
        e_total_ = totEnew
        n_accepted_scalings_ = n_accepted_scalings_ + 1
      else 
        simbox = oldbox
      end if 
    end if
  end subroutine

  !! Funktio, joka uuden ja vanhan energian perusteella
  !! p‰‰tt‰‰, hyv‰ksyt‰‰nkˆ muutos
  !!
  function acceptchange(oldenergy, newenergy) &
    result(is_accepted)
    intrinsic exp
    logical :: is_accepted
    real(dp), intent(in) :: oldenergy
    real(dp), intent(in) :: newenergy
    real(dp) :: dE
    dE = newenergy - oldenergy
    if(dE < 0._dp) then
      is_accepted = .true.
    else
      is_accepted = (grnd() < exp(-dE/temperature_))
    end if  
  end function

  !! Palauttaa hyv‰ksymisuhteet siirron ja kierron yhdistelm‰lle,
  !! sek‰ tilavuuden muutokselle
  !! 
  subroutine ratios(Nparticles, period, movratio, volratio)
    real(dp), intent(out) :: movratio, volratio
    integer, intent(in) :: Nparticles, period
    movratio = real(n_accepted_moves_, dp)/real(Nparticles*period, dp)
    volratio = real(n_accepted_scalings_, dp)/real(period, dp)   
  end subroutine

  !! P‰ivitt‰‰ maksimimuutosarvot koordinaateille/kulmille ja 
  !! sylinterin s‰teelle.
  !!
  subroutine updatemaxvalues(Nparticles, period)
    real(dp) :: mratio, vratio
    real(dp) :: olddximax, olddthetamax
    real(dp) :: newdximax, newdthetamax
    integer, intent(in) :: Nparticles, period
    !real(dp), save :: last_mratio = 0._dp
    !real(dp), save :: d_macceptance = 0._dp
    !real(dp) :: d_transrot
    call ratios(Nparticles, period, mratio, vratio)      
    !! Jos tilavuuden muutoksista on hyv‰ksytty yli 25%,
    !! kasvatetaan s‰teen maksimimuutosarvoa. Vastaavasti
    !! siirrolle/kierrolle rajana 33%.
    !! Nollataan hyv‰ksynt‰laskurit
    n_accepted_scalings_ = 0
    n_accepted_moves_ = 0
    max_scaling_ = newmaxvalue(vratio > scaling_ratio_, max_scaling_)
    call getmaxmoves(olddximax, olddthetamax)
    !if(adjust_type_ == 1 .or. adjust_type_ == 2) then
      newdthetamax = newmaxvalue(mratio > move_ratio_, olddthetamax)
      newdximax = newmaxvalue(mratio > move_ratio_, olddximax)
    !end if
    !if (adjust_type_ == 2) then
      !! adjust ratio maxtranslation/maxrotation
    !  d_macceptance = mratio - last_mratio
    !  if (d_macceptance*d_transrot > 0._dp) then
    !    newdximax = 1.05_dp*newdximax
    !    d_transrot = 1._dp
    !  else
    !    newdximax = newdximax/1.05_dp
    !    d_transrot = -1._dp
    !  end if
    !  last_mratio = mratio
    !end if
    !! :TODO: write to log if trying to increase max move too much
    !! :TODO: add also smallest possible newdximax? 
    !! :TODO: write to log if trying to decrease below smallest newdximax
    if (newdximax > largest_translation_) then
      !write(*, *) 'Tried to increase maxtrans above largest_translation.'
      newdximax = largest_translation_
      newdthetamax = largest_rotation_
    else if(newdximax < smallest_translation_) then
      !write(*, *) 'Tried to decrease maxtrans below smallest_translation.'
      newdximax = smallest_translation_
      newdthetamax = smallest_rotation_
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

  function pressure()
    real(dp) :: pressure
    pressure = pressure_
  end function

  subroutine scale_positions(oldbox, newbox, particles, n_particles)
    implicit none
    type(poly_box), intent(in) :: oldbox
    type(poly_box), intent(in) :: newbox
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(in) :: n_particles
    integer :: i
    do i = 1, n_particles
      particles(i)%x = particles(i)%x * get_x(newbox) / get_x(oldbox)
      particles(i)%y = particles(i)%y * get_y(newbox) / get_y(oldbox)
      particles(i)%z = particles(i)%z * get_z(newbox) / get_z(oldbox)
    end do    
  end subroutine

end module mc_sweep
