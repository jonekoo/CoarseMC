!Partikkelidatan ja siihen kohdistuvien operaatioiden
!määrittelyt

module particle
  use nrtype
  use gayberne, gb_save_state => save_state, gb_load_state => load_state
  use mtmod, only: grnd
  implicit none
  
  type particledat
     real(dp) :: x, y, z, ux, uy, uz
     logical :: rod
  end type particledat

  real(dp), private, save :: dthetamax, maxdr
  namelist /particle_nml/ dthetamax, maxdr
  private :: particle_nml



  contains

  subroutine save_state(write_unit)
    implicit none
    integer, intent(in) :: write_unit
    write(write_unit, NML = particle_nml)
  end subroutine save_state
  


  subroutine load_state(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    read(read_unit, NML = particle_nml)
  end subroutine load_state



  function new_particle()
    type(particledat) :: new_particle
    new_particle%x = 0._dp
    new_particle%y = 0._dp
    new_particle%z = 0._dp
    new_particle%ux = 0._dp
    new_particle%uy = 0._dp
    new_particle%uz = 1._dp
    new_particle%rod = .true.
  end function new_particle



  subroutine initParticle(maxTranslation, maxRotation)
    implicit none
    real(dp), intent(in) :: maxTranslation 
    real(dp), intent(in) :: maxRotation
    dthetamax = maxRotation
    maxdr = maxTranslation
  end subroutine initParticle



  subroutine write_module_state(unit)
    implicit none    
    integer, intent(in) :: unit
    write(unit, *) maxdr, dthetamax
  end subroutine write_module_state



  !! :TODO: put this somewhere else. Problem is that you need to initialize the potential module to use this. 
  !! 
  subroutine pairV(particlei, particlej, potE, overlap)
    implicit none
    type(particledat), intent(in) :: particlei, particlej
    real(dp), intent(out) :: potE
    logical, intent(out) :: overlap
    real(dp),dimension(3) :: rij, ui, uj
    potE = 0.0_dp
    !if(.not. (particlei%rod .and. particlej%rod)) return;
    call differences(particlei, particlej, rij(1), rij(2), rij(3))
    ui(1)=particlei%ux
    ui(2)=particlei%uy
    ui(3)=particlei%uz
    uj(1)=particlej%ux
    uj(2)=particlej%uy
    uj(3)=particlej%uz
    potE = potential(ui, uj, rij)
    overlap = .false. 
  end subroutine pairV


  
  subroutine move(oldp,newp)
    implicit none
    type(particledat), intent(in) :: oldp
    type(particledat), intent(out) :: newp
    newp = oldp
    call transmove(oldp%x,oldp%y,oldp%z,newp%x,newp%y,newp%z)
    if (oldp%rod) then
      call rotate(oldp%ux,oldp%uy,oldp%uz,newp%ux,newp%uy,newp%uz)
    end if
  end subroutine move
  

 
  subroutine transmove(xo,yo,zo,xn,yn,zn)
    use cylinder, only : getHeight
    implicit none
    real(dp), intent(in) :: xo,yo,zo
    real(dp), intent(out) :: xn,yn,zn
    real(dp) :: Lz
    xn = xo + (2.0_dp*grnd()-1.0_dp)*maxdr;
    yn = yo + (2.0_dp*grnd()-1.0_dp)*maxdr;
    zn = zo + (2.0_dp*grnd()-1.0_dp)*maxdr;
    Lz=getHeight()
    zn = zn - anint(zn/Lz)*Lz;
  end subroutine transmove



    pure subroutine differences(particlei,particlej,dx,dy,dz)
      use cylinder
      implicit none
      type(particledat),intent(in) :: particlei, particlej
      real(dp), intent(out):: dx, dy, dz
      real(dp) :: Lz
      !Juho Lintuvuori:
      !Kahden partikkelin koordinaattien erotukset
      dx = particlei%x - particlej%x;
      dy = particlei%y - particlej%y;
      dz = particlei%z - particlej%z;
      Lz = getHeight()
      !Perioidiset reunaehdot
      dz = dz - anint(dz/Lz)*Lz;
    end subroutine differences



    !Funktio, joka laskee kahden partikkelin etäisyyden
    !toisistaan
    pure function rij(particlei,particlej)
      implicit none
      type(particledat), intent(in) :: particlei,particlej
      real(dp) :: rij
      real(dp) :: dx,dy,dz,rijsq
      call differences(particlei,particlej,dx,dy,dz)      
      rijsq=dx*dx+dy*dy+dz*dz;
      rij=sqrt(rijsq) 
    end function rij



    !Funktio joka laskee partikkelien orientaatiovektorien
    !pistetulon
    real(dp) pure function idotj(particlei, particlej)
      implicit none
      type(particledat), pointer :: particlei,particlej
      idotj = particlei%ux*particlej%ux + particlei%uy*particlej%uy + &
            & particlei%uz*particlej%uz
    end function idotj
  


  !yksikkövektorin välisen pistetulon
  subroutine idotsjdots(particlei, particlej, idots, jdots)
    implicit none
    type(particledat),intent(in) :: particlei,particlej
    real(dp) :: dx,dy,dz,r,sx,sy,sz,idots,jdots
    call differences(particlei,particlej,dx,dy,dz)
    r = sqrt(dx*dx+dy*dy+dz*dz)
    sx = dx/r;
    sy = dy/r;
    sz = dz/r;
    idots = particlei%ux*sx + particlei%uy*sy + particlei%uz*sz;
    jdots = particlej%ux*sx + particlej%uy*sy + particlej%uz*sz;
  end subroutine idotsjdots



    !Palauttaa partikkelin orientaatiovektorin komponentit
    !sylinterikoordinaatistossa. 
    subroutine unitvec(particle, uro, utheta, uz)
      use utils
      implicit none
      intrinsic atan2
      type(particledat), intent(in) :: particle
      real(dp), intent(out) :: uro,utheta,uz
      real(dp) :: nx,ny,nz,uxn,uyn,uzn,theta
      theta = -atan2(particle%y, particle%x)
      nx=0.0_dp
      ny=0.0_dp
      nz=1.0_dp
      call xvec2(particle%ux, particle%uy, particle%uz, nx, ny, nz, &
               & theta, uxn, uyn, uzn)
      uro=uxn
      utheta=uyn
      uz=uzn
    end subroutine unitvec


    
  !Kierto:
  subroutine rotate(uxo, uyo, uzo, uxn, uyn, uzn)
    use utils, only : xvec2
    implicit none
    real(dp) :: uxo,uyo,uzo
    real(dp) :: uxn,uyn,uzn
    real(dp) :: theta,nx,ny,nz
    call nvec(nx,ny,nz);
    call rotangle(dthetamax,theta);
    call XVEC2(uxo,uyo,uzo,nx,ny,nz,theta,uxn,uyn,uzn)
  end subroutine rotate    



  subroutine rotangle(dthetamax,theta)
    implicit none
    double precision,intent(in) :: dthetamax
    double precision,intent(out) :: theta
    theta=(2.0_dp*grnd()-1.0_dp)*dthetamax;
  end subroutine rotangle



  subroutine nvec(nx,ny,nz)
    ! Aliohjelma satunnaisen yksikkÃ¶vektorin muodostamista varten
    ! Understanding Mol. Sim. 2nd Ed.  Frenkel, Smit
    ! s.  578
    implicit none
    double precision, intent(out) :: nx, ny, nz
    double precision :: l, u1, u2, s
    l=0.0_dp;
    do
       u1 = 1.0_dp-2.0_dp*grnd();
       u2 = 1.0_dp-2.0_dp*grnd();
       l = u1*u1+u2*u2;
       if(l <= 1.0_dp) exit;
    end do
    s = 2.0_dp*sqrt(1.0_dp-l);
    nx = u1*s;
    ny = u2*s;
    nz = 1.0_dp-2.0_dp*l;
  end subroutine nvec



  subroutine setmaxmoves(distance, angle)
    implicit none
    real(dp) :: distance,angle
    maxdr = distance
    dthetamax = angle
  end subroutine setmaxmoves



  subroutine getMaxMoves(distance, angle)
    implicit none
    real(dp), intent(out) :: distance,angle
    distance = maxdr
    angle = dthetamax
  end subroutine getMaxMoves



  !! Makes a t-configuration of particles.
  !! Comment: Let's not make box here because we may want different kind of
  !! boxes with the t-configuration inside. 
  !! 
  subroutine make_tconf(particles, n_particles)
  implicit none
  type(particledat), dimension(:), intent(inout) :: particles
  integer, intent(inout) :: n_particles
    !! Put one particle on the center of the box, parallel to z-axis
    particles(1) = new_particle()
    particles(2) = new_particle()
    !! Put other particle above that so that they form a t-configuration
    particles(2)%z = kappa_sigma()/2.0_dp + sigma_0()/2.0_dp
    particles(2)%ux = 1.0_dp
    particles(2)%uz = 0.0_dp
    n_particles = 2
  end subroutine
  


  subroutine make_sidebyside(particles, n_particles)
  implicit none
  type(particledat), dimension(:), intent(inout) :: particles
  integer, intent(inout) :: n_particles
    particles(1) = new_particle()
    particles(2) = new_particle()
    particles(1)%x = -0.5_dp*sigma_0()
    particles(2)%x = 0.5_dp*sigma_0()
    n_particles = 2
  end subroutine make_sidebyside
  
end module particle


