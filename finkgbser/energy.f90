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
  real(dp), private, save :: B0abs=11.74   !!Magnetic flux strength in Teslas
  real(dp), dimension(3), private, save :: B0vect=(/0.0,0.0,1.0/)  
  logical,private,save :: magneton=.true.
 
  contains

function diamagnetic(particlei)
  use particle
  use nrtype
  implicit none
  intrinsic dot_product
  type(particledat), intent(in) :: particlei
  real(dp) :: diamagnetic
  real(dp), dimension(3) :: B0,u,Bmol
  real(dp), dimension(3,3) :: sus
  real(dp) :: susiso,anisosus,antipar,par
  real(sp), parameter :: kB=1.38065e-23
  real(sp), parameter :: temp=107.0
  if(.not.particlei%rod) return;
  B0=B0abs*B0vect
  susiso=-3290e-30
  anisosus=1750e-30
  antipar=susiso-anisosus/3.0_dp
  par=susiso+2.0_dp*anisosus/3.0_dp
  sus(1:3,1:3)=0
  sus(1,1)=antipar
  sus(2,2)=antipar
  sus(3,3)=par
  u=(/particlei%ux,particlei%uy,particlei%uz/)
  Bmol(3)=dot_product(B0,u)
  Bmol(1)=sqrt(0.5_dp*(dot_product(B0,B0)-Bmol(3)**2))
  Bmol(2)=Bmol(1)
  diamagnetic=-dot_product(Bmol,matmul(sus,Bmol))
  diamagnetic=diamagnetic/(kB*temp)
end function



    !Palauttaa kokonaisenergian
    subroutine totenergy(prtclarray,T,ovrlp,Etot)
      implicit none
      type(particledat), dimension(:),pointer :: prtclarray
      real(dp), intent(in) :: T
      real(dp), intent(out) :: Etot
      logical, pointer :: ovrlp
      real(dp) :: Vpairtot=0.0,Vwalltot=0.0
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
      if(magneton) Etot=Etot+totDiamagnetic(prtclarray)
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
    

     
    function totDiamagnetic(particlearray)
    implicit none
    type(particledat), dimension(:),intent(in) :: particlearray
    real(dp) :: totDiamagnetic
    integer :: N,i
      totDiamagnetic=0.0
      N=size(particlearray)
      do i=1,N
        totDiamagnetic=totDiamagnetic+diamagnetic(particlearray(i))
      end do
    end function totDiamagnetic


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
      if(magneton) singleV=singleV+diamagnetic(particlei)
    end subroutine singleparticleV






  !Palauttaa hiukkasten ja seinän välisen vuorovaikutuksen
  !kokonaisenergian 
  subroutine totwallprtclV(particlearray,Eptwlltot,ovrlp)
    implicit none
    integer :: i,N,astat
    type(particledat), dimension(:),pointer :: particlearray
    real(dp), intent(out) :: Eptwlltot
    logical, pointer :: ovrlp 
    real(dp) :: oneprtclV=0.0,add=0.0,all=0.0
    type(particledat),pointer :: particlei
    real(dp), dimension(:),allocatable :: energies
    
    Eptwlltot=0.0
    ovrlp=.false.
    N=size(particlearray)
    allocate(energies(N),stat=astat)
    if(astat/=0) then
      write(*,*) 'totwallprtclV:energies'
      stop;
    end if
    energies=0.0

    do i=1,N
      particlei=>particlearray(i)
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
  end subroutine singleprtcltotV



end module energy
