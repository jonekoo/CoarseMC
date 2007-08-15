program gbcylinder
use nrtype
use mcstep, only : step,updatemaxvalues,initmcstep
use energy, only : totenergy,totwallprtclV,totpairV
use cylinder, only : initcylinder,getradius,getHeight
use verlet, only : initvlist, freevlist
use particlewall, only : initptwall,rArB
use gbgb, only : initgbgb
use particle, only : setmaxmoves
use io 
use mt19937
implicit none
! 1. muuttujien esittely
integer :: anchor      !Ankkuroinnin tyyppi 
integer,parameter :: eunit=10
integer :: avaus
character(LEN=*), parameter :: efile='energy.txt'
character(len=50) :: statefile 
type(particledat), dimension(:), pointer :: array0,ptrtoarray
type(particledat), dimension(:), allocatable,target :: particlearray
integer :: N,astat,Nrelax,Nprod,Nratio,seed,i,j,vtype
real(dp) :: radius,height,Lz,Kw,rA,rB
real(dp) :: T,pres,epses,eps0,rsphere,spmyy,epsphere,sigma0,siges,B0,B0angle
real(dp) :: totE=0.0, maxangle=0.170, maxtrans=0.156,pairE=0.0,wallE=0.0  
logical,target :: ol=.false.  
logical,pointer :: overlap
type(particledat), pointer :: particlej
!Vielä turhia muuttujia.
integer :: Nsphere
  overlap=>ol 
  !! 2. parametrien lataus
  call ReadParams(statefile,Nrelax,Nprod,Nratio,T,pres,anchor,vtype,Kw,seed, &
                  epses,eps0,rsphere,spmyy,epsphere,sigma0,siges,B0,B0angle) 
  !! 3. modulien alustus
  call initptwall(anchor,Kw)
  call init_genrand(seed)  
  call initgbgb()
  call setmaxmoves(maxtrans,maxangle)
  open(eunit,FILE=efile,status='replace',position='append',iostat=avaus)
  if(avaus/=0) then
    write (*,*) 'virhe tiedoston ',efile,' avaamisessa.'
    write (*,*) 'Ohjelman suoritus keskeytyy.'
    stop;
  end if
  !! Luetaan hiukkasten tilat tiedostosta
  call readstate(statefile,array0,radius,Lz)
  !! Alustetaan sylinteri
  call initcylinder(radius,Lz,T,pres)
  write(eunit,*) 'Sylinterin säde alussa:',getRadius()
  write(eunit,*) 'Sylinterin korkeus alussa:',getHeight()
  !! Alustetaan mcstep-moduli
  call initmcstep(vtype)
  N=size(array0)
  write (eunit,*) 'Partikkelien lukumäärä sylinterissä:', N
  allocate(particlearray(N),stat=astat)
  if(astat/=0) then
    write(*,*) 'virhe varattaessa muistia hiukkastaulukolle'
    write(*,*) 'Ohjelman suoritus keskeytyy'
    stop;
  end if
  if (associated(array0) .and. size(array0)>1) then  
     ptrtoarray=>array0
  else
    write(*,*) 'Virhe kopioitaessa partikkelitaulukkoa' 
    stop; 
  end if
  !! 6. Alustetaan Verlet'n lista (ei ehkä tässä vaan ennemmin
  !!    energiamodulissa ja tallennetaan siellä. 
  call initvlist(ptrtoarray)
  !! 7. Monte-Carlo -osuus:
  !! 7.1 Lasketaan alkutilan kokonaisenergia (overlap) (Tarvitaanko?)
  write(*,*) 'Lasketaan alkutilan energia, lämpötila=',T
  call totenergy(ptrtoarray,T,overlap,totE)
  if(overlap) then 
    write (*,*) 'gbcylinder: Alkutilassa päällekäisyyksiä.' 
    write (*,*) 'Ohjelma keskeytyy.'
    stop;
  end if
  write (eunit,*) 'Alkutilan energia: ',totE 
  !! 7.2 Relaksaatioajo (maksimimuutosarvojen säätö)
  write(*,*) 'Aloitetaan relaksaatioajo. Nrelax=',Nrelax
  do i=1,Nrelax
    call step(ptrtoarray) 
    if (mod(i,Nratio)==0) then
      call updatemaxvalues(N,Nratio)
    end if
  end do

  !! 7.3 Tuotantoajo 
  !! 7.3.1 Kirjoitetaan tiloja levylle 
  write (*,*) 'Aloitetaan tuotantoajo. Nprod=',Nprod
  do i=1,Nprod
    call step(ptrtoarray)
    if (mod(i,Nratio)==0) then
      radius=getRadius()
      height=getHeight()
      write (eunit,*) 'MC-askel:',i
      call TotEnergy(ptrtoarray,T,overlap,totE)
!      call totpairV(ptrtoarray,pairE,overlap)
!      call totwallprtclV(ptrtoarray,wallE,overlap)
!      write(eunit,*) 'Parivuorovaikutus:',pairE
!      write(eunit,*) 'GB-seinä:',wallE
      write (eunit,*) 'Kokonaisenergia:',totE
!      write (eunit,*) 'Sylinterin säde:',radius
      write (eunit,*) 'Sylinterin korkeus:',height
      call writestate(T,pres,radius,height,ptrtoarray,(i/Nratio))
    end if
  end do
   
  !! Kirjoitetaan hiukkasten loppupaikat muistiin. 
  call povout(ptrtoarray,radius,height)
  write (eunit,*) 'Ohjelman suoritus päättyi.'
  close(eunit)
  call freevlist()
  if(associated(ptrtoarray)) deallocate(ptrtoarray)
end program gbcylinder
