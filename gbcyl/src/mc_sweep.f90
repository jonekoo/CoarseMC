module mc_sweep
  use nrtype
  use energy
  use class_poly_box
  use particle
  use mtmod
  use verlet
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
  real(dp), save :: smallest_translation_ = 0.1_dp !! 1/10 molecule diameter
  real(dp), save :: largest_rotation_ = 1.57_dp     !! ~ pi/4
  real(dp), save :: smallest_rotation_ = 0.157_dp  !! 1/10 largest rotation
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
    call get_parameter(reader, 'temperature', temperature_)
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

  subroutine sweep(simbox, particles)
    type(particledat), dimension(:), intent(inout) :: particles
    type(poly_box), intent(inout) :: simbox
    integer :: i
    logical :: overlap
    integer :: i_vol_move
    real(dp) :: e_total_old
    !! Trial moves of particles and one particle move replaced with volume move
    i_vol_move = int(grnd() * real(size(particles), dp)) + 1
    e_total_old = e_total_    
    call total_energy(simbox, particles, e_total_, overlap)
    if(overlap) then 
      stop 'Overlap when entering sweep! Stopping.'
    end if
    do i = 1, size(particles)
      if (i == i_vol_move) then
        !! Replace random particle update with a volume update.
        call move_vol(simbox, particles, size(particles))
      else
        call move_particle(simbox, particles, i)
      end if
      if (mod(i, pt_ratio_) == 0) then
        call make_pt_move(simbox, particles)
      end if 
    end do
  end subroutine

  subroutine make_pt_move(simbox, particles)
    type(poly_box), intent(inout) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    real(dp) :: beta
    real(dp) :: enthalpy
    integer :: n_particles
    beta = 1._dp/temperature_
    enthalpy = e_total_ + pressure_ * volume(simbox)
    n_particles = size(particles)
    call pt_move(beta, enthalpy, particles, n_particles, simbox)
    e_total_ = enthalpy - pressure_ * volume(simbox)
    !! One way to get rid of explicit list update would be to use some kind 
    !! of observing system between the particlearray and the neighbour list. 
    call updatelist(particles, size(particles), simbox) 
  end subroutine

  subroutine move_particle(simbox, particles, i)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(in) :: i
    type(particledat) :: newparticle
    type(particledat) :: oldparticle
    logical :: overlap
    real(dp) :: e_new
    real(dp) :: e_old
    logical :: is_accepted
    e_new = 0._dp
    e_old = 0._dp
    overlap = .false.
    !! :TODO: Change moving of particles for a separate object which 
    !! :TODO: gets the whole particle array.
    call move(particles(i), newparticle)
    call set_position(newparticle, min_image(simbox, &
    (/0._dp, 0._dp, 0._dp/), position(newparticle)))
    oldparticle = particles(i) 
    particles(i) = newparticle
    call potential_energy(simbox, particles, i, e_new, overlap)
    if(.not. overlap) then 
      particles(i) = oldparticle
      call potential_energy(simbox, particles, i, e_old, overlap)
      if (overlap) then
        stop 'sweep: overlap with old particle'
      end if
      is_accepted = acceptchange(e_old, e_new)       
      if(is_accepted) then
        particles(i) = newparticle
        e_total_ = e_total_ + (e_new - e_old)
        n_accepted_moves_ = n_accepted_moves_ + 1
      end if
    end if 
  end subroutine

  !subroutine move_particleinsystem(simbox, particles, i)
  !  type(poly_box), intent(in) :: simbox
  !  type(particledat), dimension(:), intent(inout) :: particles
  !  integer, intent(in) :: i
  !  type(particledat) :: newparticle
  !  type(particledat) :: oldparticle
  !  logical :: overlap
  !  real(dp) :: e_new
  !  real(dp) :: e_old
  !  logical :: is_accepted
  !  e_new = 0._dp
  !  e_old = 0._dp
  !  overlap = .false.
    !! :TODO: Change moving of particles for a separate object which 
    !! :TODO: gets the whole particle array.
    !call move_particle(particle_iterator) 
    !call potential_energy(particle_iterator, e_new, overlap) 
    !if(.not. overlap) then 
    !  call undo_move(particle_iterator)
    !  call potential_energy(particle_iterator, e_old, overlap)
    !  if (overlap) then
    !    stop 'sweep: overlap with old particle'
    !  end if
    !  is_accepted = acceptchange(e_old, e_new)       
    !  if (is_accepted) then
    !    call redo_move(particle_iterator)
    !    e_total_ = e_total_ + (e_new - e_old) !! :TODO: Remove this.
    !    n_accepted_moves_ = n_accepted_moves_ + 1
    !  end if
    !end if 
  !end subroutine

  subroutine move_vol(simbox, particles, n_particles)
    type(poly_box), intent(inout) :: simbox
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(in) :: n_particles
    logical :: overlap
    real(dp) :: V_o, V_n
    logical :: is_accepted
    real(dp) :: tote_new
    real(dp) :: boltzmann_n
    real(dp) :: boltzmann_o
    type(poly_box) :: oldbox
    real(dp), dimension(3) :: scaling
    overlap = .false.
    V_o = volume(simbox)
    oldbox = simbox
    scaling = scale(simbox, max_scaling_, grnd)
    call scale_positions(oldbox, simbox, particles, n_particles) 
    V_n = volume(simbox)
    call total_energy(simbox, particles, tote_new, overlap)
    call scale_positions(simbox, oldbox, particles, n_particles)
    if (overlap) then
      simbox = oldbox  
    else
      boltzmann_n = tote_new + pressure_ * V_n - real(n_particles, dp) * &
        temperature_ * log(V_n)  
      boltzmann_o = e_total_ + pressure_ * V_o - real(n_particles, dp) * &
        temperature_ * log(V_o)
      is_accepted = acceptchange(boltzmann_o, boltzmann_n)
      if (is_accepted) then
        !! Scale back to new configuration
        call scale_positions(oldbox, simbox, particles, n_particles)
        e_total_ = tote_new
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
