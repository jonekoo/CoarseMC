module energy
  use nrtype, only: dp
  use particle, only: particledat
  use particlewall
!  use cell, only: list
!  use cell_energy
  !use all_pairs
  !use verlet
  use class_poly_nbrlist
  use class_poly_box
  implicit none
  private
 
  public :: totalenergy
  public :: potentialenergy

  contains
 
  subroutine totalenergy(simbox, particles, nbrlist, Etot, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_nbrlist), intent(in) :: nbrlist
    real(dp), intent(out) :: Etot
    logical, intent(out) :: overlap
    real(dp) :: Vpairtot
    real(dp) :: Vwalltot
    Etot = 0._dp
    overlap = .false.
    call pairinteractions(nbrlist, simbox, particles, Vpairtot, overlap)
    if (.not. overlap) then
      call totwallprtclV(simbox, particles, Vwalltot, overlap)
    end if  
    if (.not. overlap) then
      Etot = Vpairtot + Vwalltot
    end if
  end subroutine totalenergy
          
  !! Palauttaa hiukkasten ja seinän välisen vuorovaikutuksen
  !! kokonaisenergian 
  subroutine totwallprtclV(simbox, particles, Eptwlltot, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    real(dp), intent(out) :: Eptwlltot
    logical, intent(out) :: overlap 
    integer :: i
    real(dp) :: oneprtclV 
    Eptwlltot = 0._dp
    overlap = .false.
    !! :TODO: Could this be replaced with forall or where construct?
    do i = 1, size(particles)
      call particlewall_potential(particles(i), simbox, oneprtclV, overlap)
      if (overlap) return
      Eptwlltot = Eptwlltot + oneprtclV
    end do
  end subroutine

  !! Palauttaa yhden hiukkasen kokonaisenergian, eli
  !! vuorovaikutusenergian seinän ja muiden hiukkasten 
  !! kanssa. 
  !!
  subroutine potentialenergy(simbox, particles, nbrlist, i, Vitot, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_nbrlist), intent(in) :: nbrlist
    integer, intent(in) :: i
    real(dp), intent(out) :: Vitot
    logical, intent(out) :: overlap
    real(dp) :: Vipair 
    real(dp) :: Viwall 
    Viwall = 0._dp 
    Vipair = 0._dp
    Vitot = 0._dp
    overlap = .false.
    call particlewall_potential(particles(i), simbox, Viwall, overlap)
    if (.not. overlap) then
      call pairinteractions(nbrlist, simbox, particles, i, Vipair, overlap)
    end if
    if (.not. overlap) then 
      Vitot = Vipair + Viwall
    end if
  end subroutine

end module
