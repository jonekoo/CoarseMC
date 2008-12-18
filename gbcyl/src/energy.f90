module energy
  use nrtype, only: sp, dp
  use particle, only: particledat, rij, pairV
  use particlewall, only: prtclwallV
  use verlet, only: totpairV, singleparticleV
  implicit none
 


  public :: totenergy
  public :: singleprtcltotV
  public :: init
  public :: save_state
  public :: load_state 



  private

  real(dp), parameter :: rcut = 5.5
  real(dp), save :: magnetic_field_   !! Magnetic flux strength in Teslas
  real(dp), dimension(3), save :: magnetic_field_direction_ 
  logical, save :: is_magnet_on_  
  namelist /energy_nml/ magnetic_field_, magnetic_field_direction_, &
    & is_magnet_on_
 


  contains


 
  subroutine init(magnetic_field_direction, magnetic_field, &
    & is_magnet_on)
    implicit none
    real(dp), dimension(3), intent(in) :: magnetic_field_direction
    real(dp), intent(in) :: magnetic_field
    logical, intent(in) :: is_magnet_on
    magnetic_field_ = magnetic_field
    magnetic_field_direction_ = magnetic_field_direction
    is_magnet_on_ = is_magnet_on   
  end subroutine init
  


  subroutine save_state(write_unit)
    implicit none
    integer, intent(in) :: write_unit
    write(write_unit, NML=energy_nml)
  end subroutine save_state



  subroutine load_state(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    read(read_unit, NML=energy_nml)
  end subroutine load_state

  

  function diamagnetic(particlei)
    implicit none
    intrinsic dot_product
    type(particledat), intent(in) :: particlei
    real(dp) :: diamagnetic
    real(dp), dimension(3) :: B0,u,Bmol
    real(dp), dimension(3, 3) :: sus
    real(dp) :: susiso, anisosus, antipar, par
    real(sp), parameter :: kB=1.38065e-23
    real(sp), parameter :: temp=107.0
    if(.not.particlei%rod) return;
    B0 = magnetic_field_*magnetic_field_direction_
    susiso=-3290e-30
    anisosus=1750e-30
    antipar=susiso-anisosus/3.0_dp
    par=susiso+2.0_dp*anisosus/3.0_dp
    sus(1:3,1:3)=0
    sus(1,1)=antipar
    sus(2,2) = antipar
    sus(3,3) = par
    u=(/particlei%ux,particlei%uy,particlei%uz/)
    Bmol(3) = dot_product(B0,u)
    Bmol(1) = sqrt(0.5_dp*(dot_product(B0,B0)-Bmol(3)**2))
    Bmol(2) = Bmol(1)
    diamagnetic=-dot_product(Bmol,matmul(sus,Bmol))
    diamagnetic=diamagnetic/(kB*temp) !! :TODO: Document this!
  end function



    !Palauttaa kokonaisenergian
    subroutine totenergy(particles, n_particles, ovrlp, Etot)
      implicit none
      type(particledat), dimension(:),pointer :: particles
      integer, intent(in) :: n_particles
      real(dp), intent(out) :: Etot
      logical, intent(out) :: ovrlp
      real(dp) :: Vpairtot = 0.0, Vwalltot=0.0
      Etot=0.0
      ovrlp=.false.
      call totpairV(particles, n_particles, Vpairtot, ovrlp)
      if (ovrlp) then
        return
      end if
      call totwallprtclV(particles, n_particles, Vwalltot,ovrlp)
      if (ovrlp) then 
         return;
      end if
      Etot = Vpairtot + Vwalltot
      if(is_magnet_on_) Etot=Etot+totDiamagnetic(particles)
    end subroutine totenergy
     

     
    function totDiamagnetic(particlearray)
    implicit none
    type(particledat), dimension(:),intent(in) :: particlearray
    real(dp) :: totDiamagnetic
    integer :: n_particles,i
      totDiamagnetic=0.0
      n_particles=size(particlearray)
      do i=1,n_particles
        totDiamagnetic=totDiamagnetic+diamagnetic(particlearray(i))
      end do
    end function totDiamagnetic



  !Palauttaa hiukkasten ja seinän välisen vuorovaikutuksen
  !kokonaisenergian 
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



  !Palauttaa yhden hiukkasen kokonaisenergian, eli
  !vuorovaikutusenergian seinän ja muiden hiukkasten 
  !kanssa. 
  subroutine singleprtcltotV(particles, n_particles, particlei, i, Vitot, &
                           & ovrlp)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    logical, intent(out) :: ovrlp
    real(dp), intent(out) :: Vitot
    integer, intent(in) :: i
    real(dp) :: Vipair = 0.0, Viwall = 0.0   
    type(particledat), intent(in) :: particlei
    ovrlp=.false.
    Vitot=0.0
    call prtclwallV(particlei,Viwall,ovrlp)
    if (ovrlp) then
      return
    end if
    call singleparticleV(particles, n_particles, particlei, i, Vipair, ovrlp)
    if (ovrlp) then 
      return
    end if
    Vitot = Vipair + Viwall
  end subroutine singleprtcltotV



end module energy


