!Jatkuvasti ja tasaisesti jakautuneista LJ-partikkeleista
!koostuvan seinän ja LJ-partikkelin väliseen vuorovaikutukseen
!liittyviä parametreja ja aliohjelmia.  

module particlewall
  use cylinder,only : getRadius
  use particle, only : particledat
  use nrtype
  implicit none
  
    
  !! :TODO: Give these better names and document.
  !!
  real(dp), private, save :: alphaA, alphaB, sige, Kw
  real(dp), parameter :: LJdist=1.7 
    !LJdist = GB-partikkelin vuorovaikutuskeskuksen etäisyydet
    !partikkelin keskipisteestä.  
  logical, private, save :: ua
  namelist /particlewall_nml/ alphaA, alphaB, sige, Kw, ua
  private :: particlewall_nml



  interface 
    function repwall(x,rc,sig)
    use nrtype
    use nr, only : rf,rd
    real(dp),intent(in) :: x,rc,sig
    real(dp) :: repwall
    end function
  end interface
  
  interface
    function attwall(x,rc,sig)
    use nrtype
    use nr, only : rf,rd
    real(dp), intent(in) :: x,rc,sig
    real(dp) :: attwall
    end function attwall
  end interface 
  

contains 



  subroutine save_state(write_unit)
    implicit none
    integer, intent(in) :: write_unit
    write(write_unit, NML = particlewall_nml)
  end subroutine save_state



  subroutine load_state(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    read(read_unit, NML = particlewall_nml)
  end subroutine load_state



  !Palauttaa yhden partikkelin ja seinän välisen 
  !vuorovaikutuksen energian
  subroutine prtclwallV(prtcl,ptwallE,ovrlp)
    implicit none
    real(dp), intent(out) :: ptwallE
    logical, intent(out) :: ovrlp
    type(particledat), intent(in) :: prtcl
    ovrlp=.false.
    ptwallE=0.0
    if(prtcl%rod) then 
      call gbwall(prtcl,ptwallE,ovrlp)
    !else 
    !  call xewall(prtcl,ptwallE,ovrlp)
    end if    
  end subroutine prtclwallV



  !Palauttaa GB-partikkelin ja seinän välisen
  !vuorovaikutuksen energian. rwcut=GB-partikkelin vuorovaikutuskeskuksen
  !pienin mahdollinen etäisyys sylinterin seinästä. Syy: muuten ylivuoto
  !hypergeometrisissa
  subroutine gbwall(prtcldata, Eptwall,ovrlp)
    implicit none
    type(particledat), intent(in) :: prtcldata
    real(dp), intent(out) :: Eptwall
    logical, intent(out) :: ovrlp
    real(dp),parameter :: sig=1.0
    real(dp) :: rsiteA,rsiteB
    real(dp) :: repulA,attracA,repulB,attracB,fu,Rc
    real(dp) :: rAperR,rBperR
    real(dp),parameter :: rperRcut=0.999
    real(dp),parameter :: rwcut=0.4 
 
    ovrlp=.false.
    Eptwall=0.0   
    Rc=getradius()
    call rArB(prtcldata,rsiteA,rsiteB)
    if(rsiteA>=Rc .or. rsiteB>=Rc) then
      ovrlp=.true.
      return;
    end if
    rAperR=rsiteA/Rc
    rBperR=rsiteB/Rc
    attracA=attwall(rAperR,Rc,sig)
    repulA=repwall(rAperR,Rc,sig)
    attracB=attwall(rBperR,Rc,sig)
    repulB=repwall(rBperR,Rc,sig)
    if(.not. ua) then
      fu=1.0
    else if (ua) then
      fu=angular(prtcldata)  
    end if
    Eptwall=fu*Kw*((repulA-alphaA*attracA)+(repulB-alphaB*attracB))     
  end subroutine gbwall

  !Palauttaa partikkelin ja seinän välisten
  !vuorovaikutuskeskusten etäisyyden sylinterin keskipisteestä.
  subroutine rArB(particle,rA,rB)
    implicit none
    intrinsic sqrt
    type(particledat), intent(in) :: particle
    real(dp),intent(out) :: rA,rB
    real(dp) :: xA,yA,xB,yB
   
    xA=particle%x+LJdist*particle%ux
    yA=particle%y+LJdist*particle%uy
    xB=particle%x-LJdist*particle%ux
    yB=particle%y-LJdist*particle%uy
    rA=sqrt(xA**2 + yA**2)
    rB=sqrt(xB**2 + yB**2)
  end subroutine rArB


  !Alustaa modulin arvolla i:
  !1=random planar
  !2=uniform alignment
  !3=weak homeotropic
  subroutine initptwall(i,strength)
    implicit none
    integer,intent(in) :: i
    real(dp) :: strength
    Kw=strength
    select case(i)
    case(1)
      ua=.false.
      alphaA=1.0
      alphaB=1.0
    case(2)
      ua=.true.
      alphaA=1.0
      alphaB=1.0
    case(3)
      ua=.false.
      alphaA=1.0
      alphaB=0.0
      Kw=2*Kw      
    case default 
      write (*,*) 'Virhe alustettaessa modulia particlewall.'
      write (*,*) 'initialize(i) , i/={1,2,3}'
    end select    
  end subroutine initptwall

  real(dp) function angular(particle)
    implicit none
    type(particledat),intent(in) :: particle
    angular=(particle%uz)**2
  end function angular

end module particlewall
