!Partikkelidatan ja siihen kohdistuvien operaatioiden
!määrittelyt

module particle
  use nrtype
  use mt19937
  implicit none
  
  type particledat
     real(dp) :: x,y,z,ux,uy,uz
     logical :: rod
  end type particledat

  real(dp), save :: dthetamax,maxdr

  contains

  subroutine pairV(particlei,particlej,potE,overlap)
  use gbgb
  use nrtype
  implicit none
  type(particledat), intent(in) :: particlei, particlej
  real(dp),intent(out) :: potE
  logical, pointer :: overlap
  real(dp),dimension(3) :: rij,ui,uj
    overlap=.false.
    potE=0.0
    if(.not. (particlei%rod .and. particlej%rod)) return;
    call differences(particlei,particlej,rij(1),rij(2),rij(3))
    ui(1)=particlei%ux
    ui(2)=particlei%uy
    ui(3)=particlei%uz
    uj(1)=particlej%ux
    uj(2)=particlej%uy
    uj(3)=particlej%uz
    call gbgbV(rij,ui,uj,potE,overlap)
  end subroutine pairV


  
  subroutine move(oldp,newp)
  implicit none
  type(particledat), intent(in) :: oldp
  type(particledat), intent(out) :: newp
    newp=oldp
    call transmove(oldp%x,oldp%y,oldp%z,newp%x,newp%y,newp%z)
    if (oldp%rod) then
      call rotate(oldp%ux,oldp%uy,oldp%uz,newp%ux,newp%uy,newp%uz)
    end if
  end subroutine move
  
 
  subroutine transmove(xo,yo,zo,xn,yn,zn)
    use mt19937
    use cylinder, only : getHeight
    implicit none
  
    real(dp), intent(in) :: xo,yo,zo
    real(dp), intent(out) :: xn,yn,zn
    real(dp) :: Lz
  
    xn = xo + (2.0*genrand_real1()-1.0)*maxdr;
    yn = yo + (2.0*genrand_real1()-1.0)*maxdr;
    zn = zo + (2.0*genrand_real1()-1.0)*maxdr;
    Lz=getHeight()
    !Perioidiset reunaehdot
    zn = zn - nint(zn/Lz)*Lz;
  end subroutine transmove

    !Koordinaatistoon liittyviä funktioita

    !Palauttaa kahden partikkelin keskipisteiden etäisyydet 
    ! x-,y- ja z-suunnissa
    pure subroutine differences(particlei,particlej,dx,dy,dz)
    use cylinder
      implicit none
      type(particledat),intent(in) :: particlei,particlej
      real(dp), intent(out):: dx,dy,dz
      real(dp) :: Lz
      !Juho Lintuvuori:
      !Kahden partikkelin koordinaattien erotukset
      dx=particlei%x - particlej%x;
      dy=particlei%y - particlej%y;
      dz=particlei%z - particlej%z;
      Lz=getHeight()
      !Perioidiset reunaehdot
      dz=dz-nint(dz/Lz)*Lz;
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
      idotj=particlei%ux*particlej%ux + particlei%uy*particlej%uy + particlei%uz*particlej%uz;
    end function idotj
  

  !yksikkövektorin välisen pistetulon
  subroutine idotsjdots(particlei,particlej,idots,jdots)
      implicit none
      type(particledat),intent(in) :: particlei,particlej
      real(dp) :: dx,dy,dz,r,sx,sy,sz,idots,jdots

      call differences(particlei,particlej,dx,dy,dz)
      r=sqrt(dx*dx+dy*dy+dz*dz)
      sx=dx/r;
      sy=dy/r;
      sz=dz/r;
      idots=particlei%ux*sx + particlei%uy*sy + particlei%uz*sz;
      jdots=particlej%ux*sx + particlej%uy*sy + particlej%uz*sz;
  end subroutine idotsjdots



    !Palauttaa partikkelin orientaatiovektorin komponentit
    !sylinterikoordinaatistossa. 
    subroutine unitvec(particle,uro,utheta,uz)
      use utils
      implicit none
      intrinsic atan2
      type(particledat), intent(in) :: particle
      real(dp), intent(out) :: uro,utheta,uz
      real(dp) :: nx,ny,nz,uxn,uyn,uzn,theta
      theta=-atan2(particle%y,particle%x)
      nx=0.0
      ny=0.0
      nz=1.0
      call xvec2(particle%ux,particle%uy,particle%uz,nx,ny,nz,theta,uxn,uyn,uzn)
      uro=uxn
      utheta=uyn
      uz=uzn
    end subroutine unitvec


    

  !Kierto:
  subroutine rotate(uxo,uyo,uzo,uxn,uyn,uzn)
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
    theta=(2.0*genrand_real1()-1.0)*dthetamax;
  end subroutine rotangle



  subroutine nvec(nx,ny,nz)
    ! Aliohjelma satunnaisen yksikkÃ¶vektorin muodostamista varten
    ! Understanding Mol. Sim. 2nd Ed.  Frenkel, Smit
    ! s.  578
    implicit none
    double precision,intent(out) :: nx,ny,nz
    double precision :: l,u1,u2,s
    l=0.0;
    do
       u1=1.0-2.0*genrand_real1();
       u2=1.0-2.0*genrand_real1();
       l=u1*u1+u2*u2;
       if(l<=1.0)exit;
    end do
    s=2.0*sqrt(1.0-l);
    nx=u1*s;
    ny=u2*s;
    nz=1.0-2.0*l;
  end subroutine nvec


  subroutine setmaxmoves(distance,angle)
    implicit none
    real(dp) :: distance,angle
    maxdr=distance
    dthetamax=angle
  end subroutine setmaxmoves


  subroutine getMaxMoves(distance,angle)
    implicit none
    real(dp),intent(out) :: distance,angle
    distance=maxdr
    angle=dthetamax
  end subroutine getMaxMoves

  
end module particle





