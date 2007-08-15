!Partikkelitaulukon käsittelyyn liittyviä aliohjelmia ja muuttujia

module particledata
  use nrtype
  use particle
  use cylinder, only : rdist
  implicit none

  contains

  subroutine initptarray(particlearray,newptarray,Rc)
    implicit none
    type(particledat), dimension(:), pointer :: particlearray
    type(particledat), dimension(:), pointer :: newptarray
    type(particledat), dimension(:), pointer :: newarray
    type(particledat), dimension(:), allocatable,target :: tempptarray
    type(particledat), pointer :: particlei
    integer :: N,i,j=1,newN=0,astat
    real(dp) :: ri
    real(dp), intent(in) :: Rc
    real(dp),parameter :: cut=0.5
      
    N=size(particlearray)
    write(*,*) 'Kuutiossa on ',N,'hiukkasta'

    allocate(tempptarray(N),stat=astat)
    if(astat/=0) then 
      write (*,*) 'initptarray: Virhe varattaessa muistia tilapäiselle'
      write (*,*) 'hiukkastaulukolle. Ohjelma keskeytyy.'
      stop;
    end if
 
    do i=1,N
      particlei=>particlearray(i)
      ri=rdist(particlei%x,particlei%y)
      if( Rc-ri>cut ) then
        tempptarray(j)=particlei
        j=j+1
        newN=newN+1
      end if
    end do
    write(*,*) 'Sylinteriin tulee ',newN,'hiukkasta.'
    allocate(newarray(newN),stat=astat)
    if(astat/=0) then
      write(*,*) 'initptarray: Virhe varattaessa muistia hiukkastaulukolle.'
      write(*,*) 'Ohjelma keskeytyy.'
      stop;
    end if

    newarray=tempptarray(1:newN:1)
    !Asetetaam osoitin newptarray osoittamaan tempptarray:n alkioita
    !yhdestä newN:ään
    newptarray=>newarray 
    if (allocated(tempptarray)) then
      deallocate(tempptarray)
    end if
  end subroutine initptarray


  subroutine shell(ptrtoarray,r,dr,indicesptr,Nind)
    implicit none
    type(particledat),dimension(:),pointer :: ptrtoarray
    real(dp),intent(in) :: r,dr
    integer,dimension(:),pointer :: indicesptr,indices
    integer,dimension(:),allocatable,target :: tempind
    integer :: astat,i,j,N,Nind
    type(particledat), pointer :: particlei
    real(dp) :: Rparticlei

    N=size(ptrtoarray)
    allocate(tempind(N),stat=astat)
    if(astat/=0) then
      write(*,*) 'Virhe varattaessa muistia väliaikaiselle taulukolle'
      stop;
    end if
    j=1
    Nind=0
    do i=1,N
      particlei=>ptrtoarray(i)
      Rparticlei=rdist(particlei%x,particlei%y)
      if ( (Rparticlei<=(r+0.5*dr)) .AND. (Rparticlei>=(r-0.5*dr))) then
        tempind(j)=i
        j=j+1
        Nind=Nind+1
      end if
    end do
    if(Nind>0) then
      allocate(indices(Nind),stat=astat)
      if (astat/=0) then 
        write(*,*) 'shell:indices hyyty'
        stop;
      end if
      indices=tempind(1:j:1)
      indicesptr=>indices
      if(allocated(tempind)) deallocate(tempind)
    end if
  end subroutine

end module particledata
  
