!Verlet'n listan käsittelyyn tarvittavia muuttujia ja aliohjelmia. 
module verlet
  use nrtype
  use particle
  implicit none 
  real(dp), parameter :: rlist=6.8, rcut=5.5
  real(dp), dimension(:,:), allocatable, target, save :: xyzlist
  integer, dimension(:), allocatable, target,save :: Nvlist
  integer, dimension(:,:), allocatable,target, save ::  vlist

  contains

  !Asettaa osoittimen vlistp naapurilistaan, sekä osoittimen Nvlistp
  !naapurien lukumäärän sisältävään listan
  subroutine getvlist(vlistp,Nvlistp)
    implicit none
    integer, dimension(:,:), pointer :: vlistp
    integer, dimension(:), pointer :: Nvlistp
    
    vlistp=>vlist
    Nvlistp=>Nvlist 
  end subroutine getvlist



  !Muodostaa uuden naapurilistan hiukkastaulukosta ptarray
  subroutine newlist(ptarray)
    implicit none
    type(particledat), dimension(:), pointer :: ptarray
    type(particledat), pointer :: pti,ptj
    real(dp),dimension(:,:), pointer :: xyz
    integer :: i,N,j
    xyz=>xyzlist
    N=size(ptarray)
    do i=1,N
      Nvlist(i)=0
      xyz(1,i)=ptarray(i)%x
      xyz(2,i)=ptarray(i)%y
      xyz(3,i)=ptarray(i)%z
    end do 
    !write(*,*) 'newlist: Muodostetaan indeksilista' 
    do i=1,(N-1)
      pti=>ptarray(i)
      do j=(i+1),N
        ptj=>ptarray(j)
         if(rij(pti,ptj) < rlist ) then
           Nvlist(i)=Nvlist(i)+1
           Nvlist(j)=Nvlist(j)+1
           vlist(i,Nvlist(i))=j
           vlist(j,Nvlist(j))=i
         end if
      end do
    end do
    !write (*,*) 'newlist: Muodostettiin indeksilista'
  end subroutine newlist 



  !Päivittää naapurilistan tarvittaessa
  subroutine updatelist(ptarray)
    implicit none
    type(particledat), dimension(:), pointer :: ptarray     
    integer :: N,i,astat
    real(dp) :: dmax1=0.0, dx,dy,dz,dmax2=0.0
    type(particledat), pointer :: oldxyz,particlei

    N=size(ptarray)
    allocate(oldxyz,stat=astat)
    if (astat/=0) then
      write (*,*) 'virhe varattaessa muistia muuttujalle oldxyz'
      stop;
    end if  
    dmax1=0.0
    do i=1,N
      particlei=>ptarray(i)
      oldxyz%x=xyzlist(1,i)
      oldxyz%y=xyzlist(2,i)
      oldxyz%z=xyzlist(3,i)
      oldxyz%ux=0.0
      oldxyz%uy=0.0
      oldxyz%uz=0.0
      oldxyz%rod=.true.
      call differences(oldxyz,particlei,dx,dy,dz)
      dx=abs(dx)
      dy=abs(dy)
      dz=abs(dz)
      dmax1=max(dmax1,dx,dy,dz)
      if(dmax1 > (0.29*(rlist-rcut))) then
        call newlist(ptarray)
        exit;
      end if
    end do  
    !! dmax2=2.0*dmax1*sqrt(3.0)
    !! if(dmax2 > (rlist-rcut)) then
    !!  call newlist(ptarray)
    !! end if
    deallocate(oldxyz)
  end subroutine updatelist



  !Alustaa naapurilistan
  subroutine initvlist(ptarray)
    implicit none
    type(particledat), dimension(:), pointer :: ptarray
    integer :: astat,N,Nneigh=500

    N=size(ptarray)
    if (N<Nneigh) then
      Nneigh=N
    end if
    
    !Varataan muistia 
    allocate(xyzlist(3,N),vlist(N,Nneigh),Nvlist(N),stat=astat)
   
    if (astat/=0) then
      write (*,*) 'initvlist: Virhe varattaessa muistia.'
      write (*,*) 'Ohjelman suoritus keskeytyy.'
      stop;
    end if
    !write(*,*) 'Kutsutaan aliohjelmaa newlist aliohjelmasta initvlist'
    call newlist(ptarray)

  end subroutine initvlist

  subroutine freevlist()
  implicit none
    if (allocated(xyzlist)) deallocate(xyzlist)
    if (allocated(vlist)) deallocate(vlist)
    if (allocated(Nvlist)) deallocate(Nvlist)
  end subroutine freevlist
end module verlet
