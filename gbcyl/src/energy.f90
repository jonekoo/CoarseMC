module energy
  use nrtype, only: sp, dp
  use particle, only: particledat, rij, pairV
  use particlewall, only: prtclwallV
  use verlet, only: totpairV, singleparticleV
  implicit none
 


  public :: total_energy
  public :: potential_energy



  private

  contains
 
  !! Palauttaa kokonaisenergian
  subroutine total_energy(particles, n_particles, ovrlp, Etot)
    implicit none
    type(particledat), dimension(:), pointer :: particles
    integer, intent(in) :: n_particles
    real(dp), intent(out) :: Etot
    logical, intent(out) :: ovrlp
    real(dp) :: Vpairtot = 0.0
    real(dp) :: Vwalltot = 0.0
    Etot = 0.0
    ovrlp = .false.
    call totpairV(particles, n_particles, Vpairtot, ovrlp)
    if (.not. ovrlp) call totwallprtclV(particles, n_particles, Vwalltot,ovrlp)
    if (.not. ovrlp) then
      Etot = Vpairtot + Vwalltot
    end if
    !! write(*,*) 'Etot = ', Etot
  end subroutine total_energy
     

     
  !! Palauttaa hiukkasten ja sein�n v�lisen vuorovaikutuksen
  !! kokonaisenergian 
  subroutine totwallprtclV(particles, n_particles, Eptwlltot, ovrlp)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp), intent(out) :: Eptwlltot
    logical, intent(out) :: ovrlp 
    integer :: i
    real(dp) :: oneprtclV = 0.0
    Eptwlltot = 0.0
    ovrlp = .false.
    do i = 1, n_particles
      call prtclwallV(particles(i), oneprtclV, ovrlp)
      if (ovrlp) then 
        return;
      else
        Eptwlltot = Eptwlltot + oneprtclV
      end if
    end do
  end subroutine totwallprtclV



  !! Palauttaa yhden hiukkasen kokonaisenergian, eli
  !! vuorovaikutusenergian sein�n ja muiden hiukkasten 
  !! kanssa. 
  subroutine potential_energy(particles, n_particles, particlei, i, Vitot, &
    ovrlp)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    logical, intent(out) :: ovrlp
    real(dp), intent(out) :: Vitot
    integer, intent(in) :: i
    real(dp) :: Vipair = 0.0
    real(dp) :: Viwall = 0.0   
    type(particledat), intent(in) :: particlei
    ovrlp = .false.
    Vitot = 0.0
    call prtclwallV(particlei,Viwall,ovrlp)
    if (.not. ovrlp) then
      call singleparticleV(particles, n_particles, particlei, i, Vipair, ovrlp)
    end if
    if (.not. ovrlp) then 
      Vitot = Vipair + Viwall
    end if
  end subroutine potential_energy



end module energy


