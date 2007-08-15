module io
  use nrtype
  use particle

  contains


  !Kirjoittaa tilan eli sylinterin s‰teen, korkeuden ja partikkelitaulukon
  !tiedostoon  
  subroutine writestate(T,P,R,Lz,particlearray,index)
    implicit none
    intrinsic trim,adjustl
    real(dp), intent(in) :: T,P,R,Lz
    integer, intent(in) :: index
    type(particledat),dimension(:),pointer :: particlearray
    character(len=30) :: outfile
    integer :: GB=0,Xe=0,ios,N,astat,i
    integer, parameter :: outunit=19
    integer, dimension(:), allocatable :: help
    N=size(particlearray)
    allocate(help(N),stat=astat)
    if (astat/=0) then 
      write(*,*) 'writestate: Virhe varattaessa muistia taulukolle help.'
      stop
    end if
    GB=0
    Xe=0
    do i=1,N
      if(particlearray(i)%rod) then 
        GB=GB+1
        help(i)=1
      else 
        Xe=Xe+1
        help(i)=0
      end if
    end do 
    outfile=trim(adjustl(filename(R,T,P,GB,Xe,index)) )
    outfile=trim(outfile)
    open(outunit,file=outfile,status='replace',position='append',form='formatted',iostat=ios)
    if (ios/=0) then
      write (*,*) 'writestate: tiedostoa ',outfile,'ei voitu avata.'
      stop;
    end if
    write(outunit,*) '$R:',R,'$Lz:',Lz
    write(outunit,*) '$N:',N,'$GB:',GB,'$Xe:',Xe
    write(outunit,*) '$x:'
    write(outunit,*) particlearray%x
    write(outunit,*) '$y:'
    write(outunit,*) particlearray%y
    write(outunit,*) '$z:'
    write(outunit,*) particlearray%z
    write(outunit,*) '$rod:'
    write(outunit,*) help
    write(outunit,*) '$ux:'
    write(outunit,*) particlearray%ux
    write(outunit,*) '$uy:'
    write(outunit,*) particlearray%uy
    write(outunit,*) '$uz:'
    write(outunit,*) particlearray%uz
    close(outunit) 
    deallocate(help)
  end subroutine writestate


  subroutine readstate(filename,array0,R,Lz)
  implicit none
  character(len=*), intent(in) :: filename
  type(particledat), dimension(:), pointer :: array0
  real(dp), intent(out) :: R,Lz
  integer,parameter :: readunit=17   
  integer :: ios, N,astat,i
  character(len=3) :: charvar
  integer,dimension(:),allocatable :: help
    R=0.
    Lz=0.  
    !Koittaa avata annetun tiedoston
    open(readunit,file=filename,status='old',iostat=ios)
    if (ios/=0) then
      write(*,*) 'Tiedoston ',filename,' avaaminen ep‰onnistui.'
      write(*,*) 'Ohjelman suoritus keskeytyy.'
      stop;
    end if
    !Luetaan hiukkasten lukum‰‰r‰ tiedostosta
    read(readunit,*) charvar,R,charvar,Lz
    read(readunit,*) charvar,N
    read(readunit,*) charvar
    !Varaa muistin taulukolle
    allocate(array0(N),help(N),stat=astat)
    if (astat/=0) then
      write(*,*) 'readstate:Virhe varattaessa muistia: array0,help'
      stop;
    end if
    !Luetaan  hiukkasten x-koordinaatit
    read(readunit,*) array0(1:N)%x
    read(readunit,*) charvar
    !Luetaan y-koordinaatit
    read(readunit,*) array0(1:N)%y
    read(readunit,*) charvar,array0(1:N)%z
    !Luetaan rod-tieto
    read(readunit,*) charvar
    read(readunit,*) help(1:N)
    !Luetaan orientaatiovektoreiden komponentit
    read(readunit,*) charvar
    read(readunit,*) array0(1:N)%ux
    read(readunit,*) charvar,array0(1:N)%uy
    read(readunit,*) charvar,array0(1:N)%uz
    close(readunit) 
    do i=1,N
      if(help(i)==1) then
        array0(i)%rod=.true.
      else
        array0(i)%rod=.false.
      end if
    end do
  end subroutine readstate


  subroutine readerror(filename)
  implicit none
  character(len=*) :: filename
    write(*,*) 'Tiedoston ',filename,' muoto on v‰‰r‰.'
    write(*,*) 'Tilaa ei voitu lukea. Ohjelman suoritus keskeytyy'
    stop;
  end subroutine readerror


  !Palauttaa tiedostonimen annettujen parametrien perusteella. 
  function filename(R,T,P,GB,Xe,index) result(name)
    implicit none
    intrinsic anint
    character(len=26) :: name
    real(dp),intent(in) :: T,P,R
    integer, intent(in) :: GB,Xe,index
    integer :: Tint,Pint,Rint
    character(len=2) :: Tchar,Pchar
    character(len=3) :: Rchar
    character(len=7) :: GBchar
    character(len=5) :: Xechar
    character(len=4) :: indexchar
    integer,parameter :: mykind=selected_int_kind(8)
    Tint=anint(10*T,mykind)
    Pint=anint(10*P,mykind)
    Rint=anint(R,mykind)
    write (Tchar,'(I2.2)') Tint
    write (Pchar,'(I2.2)') Pint
    write (GBchar,'(I7.7)') GB
    write (Xechar,'(I5.5)') Xe
    write (indexchar,'(I4.4)') index
    write (Rchar,'(I3.3)') Rint
    name='R'//Rchar//'T'//Tchar//'.'//trim(adjustl(indexchar))
  end function filename



  !Aliohjelma, joka tuottaa POV-Ray -syˆttˆtiedoston partikkelitaulukosta
  !Tekij‰: Jouni Karjalainen
  !R ja Lz m‰‰r‰‰v‰t kameran ja valonl‰hteen koordinaatit. 
  subroutine povout(particlearray,R,Lz)
    implicit none
    type(particledat), dimension(:),pointer :: particlearray
    real(dp), intent(in) :: R,Lz
    integer, parameter :: povunit=18
    character(len=*),parameter :: povfile='povout.pov'
    integer :: opened,N,i
    real(dp) :: x,y,z,ux,uy,uz,camdist,litedist
    character(len=*),parameter :: GBColor='Yellow'
    character(len=*),parameter :: XeColor='Red'
    character(len=*),parameter :: CylColor='Grey'
    camdist=2.0*R
    litedist=camdist
    open(povunit,FILE=povfile,status='replace',position='append',iostat=opened)
    if(opened/=0) then
      write(*,*) 'Tiedoston ',povfile,' avaaminen ep‰onnistui.'
    else
      write(povunit,*) '#include "colors.inc"'
      write(povunit,*) '#include "transforms.inc"'
      write(povunit,*) 'background {color White}'
      write(povunit,*) 'camera {' 
      write(povunit,*) 'location <',camdist,',',camdist,',',-Lz,'>'
      write(povunit,*) 'sky <0,0,-1>'
      write(povunit,*) 'look_at <0,0,0> }'
      write(povunit,*) 'light_source { <',litedist,',',litedist,',',-Lz,'> color White}'
      write(povunit,*) 'cylinder{<0,0,',-Lz/2,'>,<0,0,',Lz/2,'>,',R,'open'
      write(povunit,*) 'texture{pigment{color ',CylColor,'}}'
      write(povunit,*) 'clipped_by{'
      write(povunit,*) 'plane{x,0}}}'
      N=size(particlearray)
      do i=1,N
        x=particlearray(i)%x
        y=particlearray(i)%y
        z=particlearray(i)%z
        ux=particlearray(i)%ux
        uy=particlearray(i)%uy
        uz=particlearray(i)%uz
        write(povunit,*) 'sphere {'
        write(povunit,*) '<0,0,0>,0.5'
        if(particlearray(i)%rod) then
          write(povunit,*) 'texture {pigment{color ',GBColor,'}}'
          write(povunit,*) 'scale <4.4,1,1>'
          write(povunit,*) 'Reorient_Trans(<1,0,0>,<',ux,',',uy,',',-uz,'>)'
        else
          write(povunit,*) 'scale 0.96'
          write(povunit,*) 'texture {pigment{color ',XeColor,'}}'
        end if
        write(povunit,*) 'translate <',x,',',y,',',-z,'>}'
      end do  
      close(povunit)
    end if
  end subroutine povout



  !Aliohjelma parameterien lukemiseen tiedostosta paramsfile.
  !Tekij‰: Juho Lintuvuori
  !muutokset: Jouni Karjalainen
  subroutine ReadParams(file,Nrelax,Nprod,Nratio,T,pres,anchor,voltyp,Kw,&
                        & seed,epses,eps0,rsphere,spmyy,epsphere,sigma0,siges,&
                        & B0, B0angle, domainw, maxdr, cutoff)
    implicit none
    ! Alkuper√§inen aliohjelma Antti Kurosen k√§sialaa; kurssilta johdatus
    ! atomistisiin simulaatioihin
    integer,intent(out) :: Nrelax,Nprod,Nratio,anchor,seed
    integer,intent(out) :: voltyp
    real(dp),intent(out) :: T,pres,epses,eps0,rsphere,spmyy,epsphere 
    real(dp),intent(out) :: sigma0,siges,Kw,B0,B0angle,domainw,maxdr,cutoff
    character(len=*), parameter :: paramsfile='gbcyl.in'
    integer, parameter :: input=20
    character(len=50),intent(out) :: file  
    integer :: i, ios
    character(len=100) :: line
    character(len=50) :: string
    real(dp) :: x
    open(input,file=paramsfile,status='old',iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Open error for file ',paramsfile,'.' 
      STOP 
    end if
    ! Loop through file and find all variables
    i=0
    do
     i=i+1
     ! First read in line in whole.
     read(input,fmt='(A80)',iostat=ios) line
     ! Check for probable end-of-file
     if (ios < 0) exit
     ! If not parameter line, cycle
     if (line(1:1) /= '$') cycle
     if (line(1:10)=='$initstate') then
       read(unit=line,fmt=*,iostat=ios) string,file
     else       
       ! Attempt to parse line into command and value
       ! using internal formatted read
       read(unit=line,fmt=*,iostat=ios) string,x
       if (ios > 0) then
          print *,'ERROR: Data read in error on line',i
          stop 'Invalid parameter file'
       endif
       ! From here on line should be in variable format
       ! Look for variable-defining strings
       if (string=='$Nrelax') then
	 Nrelax=int(x+0.5)
       else if (string=='$Nprod') then
	 Nprod=int(x+0.5);
       else if (string=='$Nratio') then
	 Nratio=int(x+0.5);
       else if (string=='$T') then
	 T=x 
       else if (string=='$pres') then
	 pres=x
       else if (string=='$anchor') then
	 anchor=int(x+0.5)
       else if (string=='$voltyp') then
         voltyp=int(x+0.5)
       else if (string=='$seed') then
         seed=int(x+0.5)
       else if(string=='$epses')then
         epses=x;
       else if(string=='$eps0')then
         eps0=x;
       else if(string=='$rsphere')then
         rsphere=x;
       else if(string=='$spmyy')then
         spmyy=x;
       else if(string=='$epsphere')then
         epsphere=x;
       else if(string=='$sigma0')then
         sigma0=x;
       else if(string=='$siges')then
         siges=x;
       else if(string=='$Kw')then
         Kw=x;
       else if(string=='$B0')then
         B0=x;
       else if(string=='$B0angle')then
         B0angle=x;
       else if(string=='$cutoff')then
         cutoff=x;
       else if(string=='$domainw')then
         domainw=x;
       else if(string=='$maxdr')then
         maxdr=x;
       else
         print '(A,A)','Unknown parameter',string
	 stop 'Parameter read in error'
       endif
     end if
       print '(A,A16,A,G13.6)','Read in parameter ',string,' value',x
    
    enddo

    close(input)
  
  
  end subroutine ReadParams

  
  

end module io
