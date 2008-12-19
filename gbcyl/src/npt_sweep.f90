  subroutine sweep(particles, n_particles)
    use particle
    use cylinder
    use mc_sweep
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



