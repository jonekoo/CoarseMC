module energy
  use nrtype
  use particle
  use gbgb
  use particlewall
  use verlet,only : getvlist
  use cylinder, only : volume,getPres,getradius
  implicit none
  
  real(dp), parameter :: rcut=5.5
  real(dp), private, save :: totwallE
  real(dp), private, save :: B0abs=11.74   !! Magnetic flux density in Teslas
  real(dp), dimension(3), private, save :: B0vect=(/0.0,0.0,1.0/)  
  logical,private,save :: magneton=.false.
 
  contains

    
    subroutine setB0Field(B0, B0theta, B0phi)
    use nrtype
    implicit none
    real(dp), intent(in) :: B0, B0theta, B0phi
      magneton=.true.
      B0abs=B0
      B0vect(1)=sin(B0theta)*cos(B0phi)
      B0vect(2)=sin(B0theta)*sin(B0phi)
      B0vect(3)=cos(B0theta)
    end subroutine setB0Field



    !Palauttaa kokonaisenergian
    subroutine totenergy(prtclarray,T,ovrlp,Etot)
      implicit none
      type(particledat), dimension(:),pointer :: prtclarray
      real(dp), intent(in) :: T
      real(dp), intent(out) :: Etot
      logical, pointer :: ovrlp
      real(dp) :: Vpairtot=0.0, Vwalltot=0.0
      integer :: i,N
      Etot=0.0
      ovrlp=.false.
      !write(*,*) 'Lasketaan kokonaisparivuorovaikutus'
      call totpairV(prtclarray,Vpairtot,ovrlp)
      if (ovrlp) then
        !write (*,*) 'Parivuorovaikutuksessa päällekäisyys'
        return
      end if
      !write (*,*) 'Lasketaan kokonaisvuorovaikutus seinän kanssa'
      call totwallprtclV(prtclarray,Vwalltot,ovrlp)
      if (ovrlp) then 
        !write(*,*) 'Päällekäisyys seinän kanssa'
        return;
      end if
      Etot= Vpairtot + Vwalltot
      if(magneton) Etot=Etot+totMagnetic(prtclarray)
    end subroutine totenergy
 


    function pv(T,N) result(E)
      implicit none
      real(dp) :: V,P
      real(dp) :: E
      real(dp), intent(in) :: T
      integer, intent(in) :: N
     
      E=0.0
      V=volume()
      P=getPres()
      E=P*V-real(N)*T*log(V)
    end function pv


    !Palauttaa parivuorovaikutuksen kokonaisenergian
    subroutine totpairV(particlearray,Vtot,ovrlp)
      implicit none
      integer :: i,N,j,jj
      type(particledat), dimension(:), pointer :: particlearray
      real(dp), intent(out) :: Vtot
      real(dp) :: pairE
      logical, pointer :: ovrlp
      integer,dimension(:,:), pointer :: verletl
      integer,dimension(:), pointer :: nvlist
      type(particledat),pointer :: particlei,particlej
     
      N=size(particlearray)      
      Vtot=0.0
      call getvlist(verletl,nvlist)
      do i=1,N
        if(nvlist(i)==0) cycle;
        particlei=>particlearray(i)
        do jj=1,nvlist(i)
          j=verletl(i,jj)
          if(j<=i) cycle;
          particlej=>particlearray(j)
          if (rij(particlei,particlej)<rcut) then 
            call pairV(particlei,particlej,pairE,ovrlp)
            if (ovrlp) return;
            Vtot=Vtot+pairE      
          else 
            cycle
          end if 
        end do
      end do
    end subroutine totpairV
    

     
    function totMagnetic(parray)
    use particle, only : particledat, magnetic
    implicit none
    type(particledat), dimension(:),intent(in) :: parray
    real(dp) :: totMagnetic
    integer :: N,i
      totMagnetic=0.0
      N=size(parray)
      do i=1,N
        totMagnetic=totMagnetic+magnetic(B0abs,B0vect,parray(i))
      end do
    end function totMagnetic


    !Laskee hiukkasen i vuorovaikutusenergian muiden hiukkasten
    !kanssa
    subroutine singleparticleV(prtclarray,particlei,i,singleV,ovrlp)
    implicit none
    integer :: i,j,N
    type(particledat),dimension(:),pointer :: prtclarray
    integer, dimension(:,:), pointer :: verletl
    integer, dimension(:), pointer :: Nverletl
    real(dp) :: singleV,pairE=0.0
    logical, pointer :: ovrlp
    type(particledat) :: particlei,particlevij    
      call getvlist(verletl,Nverletl)
      N=Nverletl(i)
      singleV=0.0
      do j=1,N
        particlevij=prtclarray(verletl(i,j))
        if (rij(particlei,particlevij)<rcut) then 
          call pairV(particlei,particlevij,pairE,ovrlp)
          if (ovrlp) return;
          singleV=singleV+pairE
        else 
          cycle;
        end if
      end do
    end subroutine singleparticleV






  !Palauttaa hiukkasten ja seinän välisen vuorovaikutuksen
  !kokonaisenergian 
  subroutine totwallprtclV(parray,Eptwlltot,ovrlp)
    implicit none
    integer :: i,N,astat
    type(particledat), dimension(:),pointer :: parray
    real(dp), intent(out) :: Eptwlltot
    logical, pointer :: ovrlp 
    real(dp) :: oneprtclV=0.0,add=0.0,all=0.0
    type(particledat),pointer :: particlei
    real(dp), dimension(:),allocatable :: energies
    
    Eptwlltot=0.0
    ovrlp=.false.
    N=size(parray)
    allocate(energies(N),stat=astat)
    if(astat/=0) then
      write(*,*) 'totwallprtclV:energies'
      stop;
    end if
    energies=0.0

    do i=1,N
      particlei=>parray(i)
      call prtclwallV(particlei,oneprtclV,ovrlp)
      if (ovrlp) then 
        return;
      else
        energies(i)=oneprtclV
      end if
    end do
    Eptwlltot=sum(energies)
    deallocate(energies)

  end subroutine totwallprtclV


  !Palauttaa yhden hiukkasen kokonaisenergian, eli
  !vuorovaikutusenergian seinän ja muiden hiukkasten 
  !kanssa. 
  subroutine singleprtcltotV(prtclarray,particlei,i,Vitot,ovrlp)
    implicit none
    logical, pointer :: ovrlp
    real(dp), intent(out) :: Vitot
    type(particledat), dimension(:), pointer :: prtclarray
    integer, intent(in) :: i
    real(dp) :: Vipair=0.0,Viwall=0.0   
    type(particledat), intent(in) :: particlei
    ovrlp=.false.
    Vitot=0.0
    call prtclwallV(particlei,Viwall,ovrlp)
    if (ovrlp) then
      return
    end if
    call singleparticleV(prtclarray,particlei,i,Vipair,ovrlp)
    if (ovrlp) then 
      return
    end if
    Vitot=Vipair+Viwall
    if(magneton) Vitot=Vitot+magnetic(B0abs,B0vect,particlei)
  end subroutine singleprtcltotV



end module energy
