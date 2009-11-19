module energy
  use nrtype, only: dp
  use particle, only: particledat
  use particlewall
  use verlet
  use class_poly_box
  implicit none
  private
 
  public :: total_energy
  public :: potential_energy

  contains
 
  !! Palauttaa kokonaisenergian
  subroutine total_energy(particles, n_particles, simbox, Etot, overlap)
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: Etot
    logical, intent(out) :: overlap
    real(dp) :: Vpairtot
    real(dp) :: Vwalltot
    Etot = 0._dp
    overlap = .false.
    call pair_interactions(particles, n_particles, simbox, Vpairtot, overlap)
    if (.not. overlap) then
      call totwallprtclV(particles, n_particles, simbox, Vwalltot, overlap)
    end if  
    if (.not. overlap) then
      Etot = Vpairtot + Vwalltot
    end if
  end subroutine total_energy
          
  !! Palauttaa hiukkasten ja seinän välisen vuorovaikutuksen
  !! kokonaisenergian 
  subroutine totwallprtclV(particles, n_particles, simbox, Eptwlltot, overlap)
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: Eptwlltot
    logical, intent(out) :: overlap 
    integer :: i
    real(dp) :: oneprtclV 
    Eptwlltot = 0._dp
    overlap = .false.
    do i = 1, n_particles
      call particlewall_potential(particles(i), simbox, oneprtclV, overlap)
      if (overlap) then 
        return;
      else
        Eptwlltot = Eptwlltot + oneprtclV
      end if
    end do
  end subroutine totwallprtclV

  !! Palauttaa yhden hiukkasen kokonaisenergian, eli
  !! vuorovaikutusenergian seinän ja muiden hiukkasten 
  !! kanssa. 
  subroutine potential_energy(particles, n_particles, particlei, i, simbox, &
    Vitot, overlap)
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    type(particledat), intent(in) :: particlei
    integer, intent(in) :: i
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: Vitot
    logical, intent(out) :: overlap
    real(dp) :: Vipair = 0._dp
    real(dp) :: Viwall = 0._dp 
    Vitot = 0._dp
    overlap = .false.
    call particlewall_potential(particlei, simbox, Viwall, overlap)
    if (.not. overlap) then
      call pair_interactions(particles, n_particles, simbox, particlei, i, &
        Vipair, overlap)
    end if
    if (.not. overlap) then 
      Vitot = Vipair + Viwall
    end if
  end subroutine potential_energy

end module energy

