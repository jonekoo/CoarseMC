module io
  use nrtype
  use particle

  
  real(dp) :: radius,height
  integer :: partccount,GBcount,Xecount
  namelist /gbxecyl/ radius, height, partccount, GBcount, Xecount   


  interface readState
    module procedure readState, readStateOld  
  end interface readState


  contains

  subroutine writeState(T,P,R,Lz,N,parray,index)
  use particle, only : particledat, particleToString 
  implicit none
  real(dp), intent(in) :: T,P,R,Lz
  integer, intent(in) :: N,index
  type(particledat), dimension(:), intent(in) :: parray
  integer :: i,ios,GB,Xe
  integer, parameter :: outunit=19
  character(len=30) :: outfile !! ="tempstatefile.out"
  character(len=140) :: string
    GB=N
    Xe=0
    outfile=trim(adjustl(filename(R,T,P,GB,Xe,index)))
    open(outunit,file=outfile,status='replace', position='append', &
       & form='formatted', iostat=ios)
    if(ios/=0) then 
      write(*,*) 'writestate: tiedostoa ',outfile,' ei voitu avata.'
      stop;
    end if
    radius=R
    height=Lz
    partccount=N
    GBcount=GB
    Xecount=Xe
    write(outunit,NML=gbxecyl)
    do i=1,N
      call particleToString(parray(i), string)
      write(outunit,*) trim(adjustl(string))
    end do
    close(outunit)
  end subroutine writeState


  subroutine readState(filename,N,parray,R,Lz)
  use particle, only : particledat, stringToParticle
  use nrtype
  implicit none
  character(len=*), intent(in) :: filename  
  integer, intent(out) :: N
  type(particledat), dimension(:), pointer :: parray
  real(dp), intent(out) :: R,Lz
  integer, parameter :: readunit=17
  integer :: ios,i
  character(len=140) :: string
    open(readunit,file=filename,status='old',iostat=ios)
    if(ios/=0) then 
      write(*,*) 'Failed to open statefile ',filename
      stop;
    end if
    read(readunit,nml=gbxecyl)
    N=partccount
    !GB=GBcount
    !Xe=Xecount
    R=radius
    Lz=height
    allocate(parray(N))
    do i=1,N
      read(readunit,'(A)',iostat=ios) string
      if(ios<0)then
        write(*,*) 'partikkeleita luettiin',i  
        stop;
      end if
      call stringToParticle(string,parray(i))
!      write(*,*) trim(adjustl(string))
    end do
    close(readunit)
  end subroutine readState



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
  subroutine ReadParams(file,filetype,Nrelax,Nprod,Nratio,T,pres,anchor,&
                     &  voltyp,Kw, seed,epses,eps0,rsphere,spmyy,epsphere,&
                     &  sigma0,siges, B0, B0angle, magneton, domainw, maxdr,&
                     &  cutoff)
    implicit none
    ! Alkuper√§inen aliohjelma Antti Kurosen k√§sialaa; kurssilta johdatus
    ! atomistisiin simulaatioihin
    integer,intent(out) :: Nrelax,Nprod,Nratio,anchor,seed,filetype
    integer,intent(out) :: voltyp, magneton
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
       else if(string=='$magneton')then
         magneton=int(x+0.5)
       else if(string=='$cutoff')then
         cutoff=x;
       else if(string=='$domainw')then
         domainw=x;
       else if(string=='$maxdr')then
         maxdr=x;
       else if(string=='$statefiletype')then
         filetype=int(x+0.5);
       else
         print '(A,A)','Unknown parameter',string
	 stop 'Parameter read in error'
       endif
     end if
       print '(A,A16,A,G13.6)','Read in parameter ',string,' value',x
    
    enddo
    close(input)
  end subroutine ReadParams


  
  subroutine readStateOld(filename,N,xs,ys,zs,uxs,uys,uzs,rods,R,Lz)
  use nrtype
  implicit none
  character(len=*), intent(in) :: filename
  integer,intent(out) :: N
  real(dp),dimension(:),pointer :: xs,ys,zs,uxs,uys,uzs
  logical,dimension(:),pointer :: rods
  real(dp), intent(out) :: R,Lz
  integer,parameter :: readunit=17   
  integer :: ios,astat,i
  character(len=5) :: charvar
  integer,dimension(:),allocatable :: help
  integer :: GB,Xe
  character(len=30) :: errorstring="statefile format is incompatible"
  real(dp) :: radius,height
  integer :: partccount,GBcount,Xecount
  namelist /gbxecyl/ radius, height, partccount, GBcount, Xecount   
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
    read(readunit,nml=gbxecyl) 
    N=partccount;  !write(*,*) N
    GB=GBcount;    !write(*,*) GB
    Xe=Xecount;    !write(*,*) Xe
    R=radius;      !write(*,*) R
    Lz=height;     !write(*,*) Lz
    if (R<=0) then 
      write(*,*) 'readstate: Cylinder radius is not a positive number.'
      stop;
    end if
    if(N/=(GB+Xe)) then
      write(*,*) 'readstate: GB+Xe/=N, particlecount doesnt match in &
                & the file',filename
      stop;
    end if
    !Varaa muistin taulukolle
    allocate(xs(N),ys(N),zs(N),uxs(N),uys(N),uzs(N),rods(N),help(N),stat=astat)
    if (astat/=0) then
      write(*,*) 'readstate:Virhe varattaessa muistia'
      stop;
    end if
    !Luetaan  hiukkasten x-koordinaatit
    read(readunit,*) charvar;  !write(*,*) charvar
    read(readunit,*) xs(1:N);  
    read(readunit,*) charvar; !write(*,*) charvar
    !if(charvar(1:1)/='$') write(*,*) errorstring; stop;  
    !Luetaan y-koordinaatit
    read(readunit,*) ys(1:N)   
    read(readunit,*) charvar
    read(readunit,*) zs(1:N)
    !Luetaan rod-tieto
    read(readunit,*) charvar
    read(readunit,*) help(1:N)
    !Luetaan orientaatiovektoreiden komponentit
    read(readunit,*) charvar
    read(readunit,*) uxs(1:N)
    read(readunit,*) charvar
    read(readunit,*) uys(1:N)
    read(readunit,*) charvar
    read(readunit,*) uzs(1:N)
    close(readunit)    
    do i=1,N
      if(help(i)==1) then
        rods(i)=.true.
      else
        rods(i)=.false.
      end if
    end do
  end subroutine readStateOld
  
  
end module io
