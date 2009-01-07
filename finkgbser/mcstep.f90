module mcstep
  use nrtype
  use energy
  use cylinder
  use particle
  use mt19937, only : genrand_real1
  use verlet, only : updatelist

  integer,private,save :: accmov=0,accvol=0 
  real,private,save :: movrat,volrat
  real(dp),private,save :: maxdRadius=0.1
  integer, private, save :: voltyp 
  
  contains

    subroutine step(prtclarray)
      use cylinder, only : getTemp,getPres
      implicit none
      intrinsic log,real
      type(particledat), dimension(:),pointer :: prtclarray
      integer :: i,N
      real(dp) :: Eold=0.0, Enew=0.0, totEold=0.0, totEnew=0.0
      real(dp) :: x0,y0,z0,xn,yn,zn
      real(dp) :: ux0,uy0,uz0,uxn,uyn,uzn
      real(dp) :: R0,Rn,Lz0,Lzn
      real(dp) :: totalenergy=0.0
      real(dp) :: V0,Vn
      logical :: accept=.true.
      logical,target :: ol=.false.
      logical, pointer :: overlap
      type(particledat) :: newparticle,oldparticle 
      real(dp) :: tempr,press
      tempr=getTemp()
      press=getPres()      
      !write(*,*) 'Aloitetaan hiukkasten yksitt‰inen siirto ja kierto.'
      N=size(prtclarray)
      !write(*,*) 'Hiukkasten lkm: ',N
      overlap=>ol
      do i=1,N
        overlap=.false.
        oldparticle=prtclarray(i)
        call move(oldparticle,newparticle)
        call singleprtcltotV(prtclarray,newparticle,i,Enew,overlap)
        if(overlap) then 
          prtclarray(i)=oldparticle
          cycle;
        end if 
        call singleprtcltotV(prtclarray,oldparticle,i,Eold,overlap)   
        accept=acceptchange(Eold,Enew,tempr,overlap)       
        if(.not. accept) then     
          prtclarray(i)=oldparticle
        else
          prtclarray(i)=newparticle
          accmov=accmov+1
        end if
      end do
      !Tallennetaan vanha tilavuus
      V0=volume()
!      if (voltyp==1) then
!        !write (*,*) 'Tallennetaan kokonaisenergia ennen tilavuuden muutosta' 
!        call totenergy(prtclarray,tempr,overlap,totEold)
!        R0=getRadius()
!        Rn=R0+(2.0*genrand_real1()-1.0)*maxdRadius
!        call changeR(prtclarray,R0,Rn);
!        call totenergy(prtclarray,tempr,overlap,totEnew)
!      else if(voltyp==2) then
        !S‰de ei muutu, joten partikkeli-sein‰ -vuorovaikutus vakio
        !tilavuuden muutoksessa. 
        Lz0=getHeight()
        Lzn=Lz0+(2.0*genrand_real1()-1.0)*maxdRadius
        call changeLz(prtclarray,Lz0,Lzn);
        Vn=volume()
        call totpairV(prtclarray,totEnew,overlap)
        call changeLz(prtclarray,Lzn,Lz0);
        if(overlap) return; !! If the new state results in an overlap, 
                            !! old state is restored
        call totpairV(prtclarray,totEold,overlap)  !! Total energy for the old
                                                   !! state
!      end if
      totEnew=totEnew+press*Vn-real(N)*tempr*log(Vn)   
      totEold=totEold+press*V0-real(N)*tempr*log(V0)
      accept=acceptchange(totEold,totEnew,tempr,overlap)
      if(accept) then
!       if(voltyp==1) then 
!          call changeR(prtclarray,Rn,R0);
!        else if(voltyp==2) then
          call changeLz(prtclarray,Lz0,Lzn);
!        end if
        accvol=accvol+1
      end if 
      !write (*,*) 'P‰ivitet‰‰n Verlet-lista, jos tarpeen.'
      call updatelist(prtclarray)
    end subroutine step




    !Funktio, joka uuden ja vanhan energian perusteella
    !p‰‰tt‰‰, hyv‰ksyt‰‰nkˆ muutos
    logical function acceptchange(oldenergy,newenergy,tempr,overlap) result(ok)
      implicit none
      intrinsic exp
      real(dp),intent(in) :: oldenergy, newenergy, tempr
      real(dp) :: u, dE
      logical, intent(in) :: overlap
      
      if(overlap) then 
        ok=.false.
        return;
      end if
      dE=newenergy-oldenergy
      ok=.false.
      if(dE<0) then
        ok=.true.
      else
        u=genrand_real1()
        if ( u<exp(-dE/tempr) ) then 
          ok=.true.
        end if
      end if  
    end function acceptchange



    !Palauttaa hyv‰ksymisuhteet siirron ja kierron yhdistelm‰lle,
    !sek‰ tilavuuden muutokselle
    subroutine ratios(Nparticles,period,movratio,volratio)
      implicit none
      real,intent(out) :: movratio,volratio
      integer, intent(in) :: Nparticles,period

      movratio=real(accmov/(Nparticles*period))
      volratio=real(accvol/period)
      movrat=movratio
      volrat=volratio
    end subroutine ratios



    !P‰ivitt‰‰ maksimimuutosarvot koordinaateille/kulmille ja 
    !sylinterin s‰teelle.
    subroutine updatemaxvalues(Nparticles,period)
      implicit none
      real :: mratio,vratio
      logical :: volinc,movinc
      real(dp) :: olddximax,olddthetamax
      real(dp) :: newdRmax,newdximax,newdthetamax
      integer,intent(in) :: Nparticles,period

      call ratios(Nparticles,period,mratio,vratio)      
      !Jos tilavuuden muutoksista on hyv‰ksytty yli 25%,
      !kasvatetaan s‰teen maksimimuutosarvoa. Vastaavasti
      !siirrolle/kierrolle rajana 33%.
      volinc=(0.25<vratio)
      movinc=(0.33<mratio)
      
      !Nollataan hyv‰ksynt‰laskurit
      accvol=0
      accmov=0

      newdRmax=newmaxvalue(volinc,maxdRadius)      
      maxdRadius=newdRmax

      call getmaxmoves(olddximax,olddthetamax)
      newdthetamax=newmaxvalue(movinc,olddthetamax)
      newdximax=newmaxvalue(movinc,olddximax)
      call setmaxmoves(newdximax,newdthetamax)
  
    end subroutine updatemaxvalues
  

    function newmaxvalue(increase,oldvalue) result(newvalue)
      implicit none
      logical, intent(in) :: increase
      real(dp) :: newvalue
      real(dp),intent(in) :: oldvalue

      newvalue=oldvalue
      if(increase) then
        newvalue=oldvalue*1.05
      else 
        newvalue=oldvalue*0.95
      end if
    end function newmaxvalue
    

    subroutine changeR(particlearray, oldR,newR)
      implicit none
      real(dp),intent(in) :: oldR,newR
      type(particledat), dimension(:), pointer :: particlearray 
      integer :: i
      real(dp) :: x0,y0        

      call setR(newR)
      do i=1,size(particlearray)
        x0=particlearray(i)%x
        y0=particlearray(i)%y
        call scale(x0,y0,particlearray(i)%x,particlearray(i)%y,oldR,newR)
      end do
    end subroutine changeR


    subroutine changeLz(particlearray,oldLz,newLz)
      implicit none
      real(dp), intent(in) :: oldLz,newLz
      type(particledat),dimension(:),pointer :: particlearray
      integer :: i
      real(dp) :: z0
      
      call setLz(newLz)
      do i=1,size(particlearray)
        particlearray(i)%z=newLz/oldLz*particlearray(i)%z
      end do
    end subroutine changeLz

    subroutine initmcstep(vtype)
      implicit none
      integer, intent(in) :: vtype
      voltyp=vtype
    end subroutine initmcstep

end module mcstep
