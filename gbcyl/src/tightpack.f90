program tightpack
  use nrtype
  use particle, only : particledat
  use particledata, only : initptarray
  use io, only : writestate
  integer :: N,Nalloc=104544,astat
  type(particledat),dimension(:),allocatable,target :: particlearray
  type(particledat),dimension(:),pointer :: array0,temparray
  real(dp) :: Lx,Ly,Lz
  real(dp), parameter :: T=2.0, P=2.0,Rc=60.0
  integer,parameter :: index=000
 
  !Varataan muistia hiukkastaulukolle
  allocate(particlearray(Nalloc),stat=astat)
  if(astat/=0) then
    write(*,*) 'Tightpack: Virhe varattaessa muistia: particlearray'
    stop;
  end if
  
  !Luetaan sis‰‰n partikkelikuutio, jossa 
  !heksagonaalinen tiivispakkaus.  
  call readinatoms(N,Lx,Ly,Lz,particlearray)
 
  array0=>particlearray  

  !Valitaan hiukkasista sylinterin sis‰‰n mahtuvat.
  call initptarray(array0,temparray,Rc)
  deallocate(particlearray)  
  !Kirjoitetaan uusi taulukko tiedostoon. 
  call writestate(T,P,Rc,Lz,temparray,index)
 

end program tightpack
