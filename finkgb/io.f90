module io
  use nrtype
  implicit none

  real(dp) :: radius,height
  integer :: partccount,GBcount,Xecount
  namelist /gbxecyl/ radius, height, partccount, GBcount, Xecount   
  contains
  

  !Kirjoittaa tilan eli sylinterin s‰teen, korkeuden ja partikkelitaulukon
  !tiedostoon  
  subroutine writestate(T,P,R,Lz,N,xs,ys,zs,uxs,uys,uzs,rods,index)
    implicit none
    intrinsic trim,adjustl
    real(dp), intent(in) :: T,P,R,Lz
    integer, intent(in) :: index
    real(dp),dimension(:),intent(in) :: xs,ys,zs,uxs,uys,uzs
    logical,dimension(:),intent(in) :: rods
    integer,intent(in) :: N
    character(len=30) :: outfile
    integer :: GB=0,Xe=0,ios,astat,i
    integer, parameter :: outunit=19
    integer, dimension(:), allocatable :: help

    allocate(help(N),stat=astat)
    if (astat/=0) then 
      write(*,*) 'writestate: Virhe varattaessa muistia taulukolle help.'
      stop
    end if

    GB=0
    Xe=0
    do i=1,N
      if(rods(i)) then 
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
    radius=R
    height=Lz
    partccount=N
    GBcount=GB
    Xecount=Xe
    write(outunit,NML=gbxecyl)
    write(outunit,*) '$x:'
    write(outunit,*) xs(1:N)
    write(outunit,*) '$y:'
    write(outunit,*) ys(1:N)
    write(outunit,*) '$z:'
    write(outunit,*) zs(1:N)
    write(outunit,*) '$rod:'
    write(outunit,*) help(1:N)
    write(outunit,*) '$ux:'
    write(outunit,*) uxs(1:N)
    write(outunit,*) '$uy:'
    write(outunit,*) uys(1:N)
    write(outunit,*) '$uz:'
    write(outunit,*) uzs(1:N)
    close(outunit) 
    deallocate(help)
  end subroutine writestate


  subroutine readstate(filename,N,xs,ys,zs,uxs,uys,uzs,rods,R,Lz)
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
!    read(readunit,*) charvar,R,charvar,Lz
!    read(readunit,*) charvar,N,charvar,GB,charvar,Xe
    if (R<=0) then 
      write(*,*) 'readstate: Cylinder radius is not a positive number.'
      stop;
    end if
    if(N/=(GB+Xe)) then
      write(*,*) 'readstate: GB+Xe/=N, particlecount doesnt match in &
                 the file',filename
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






  !Aliohjelma parameterien lukemiseen tiedostosta paramsfile.
  !Tekij‰: Juho Lintuvuori
  !muutokset: Jouni Karjalainen
  subroutine ReadParams(file,Nrelax,Nprod,Nratio,T,pres,anchor,voltyp,Kw,&
                        & seed,epses,eps0,rsphere,spmyy,epsphere,sigma0,siges,&
                        & allign,debug,cutoff,maxdr,domainw)
    implicit none
    ! Alkuper√§inen aliohjelma Antti Kurosen k√§sialaa; kurssilta johdatus
    ! atomistisiin simulaatioihin
    integer,intent(out) :: Nrelax,Nprod,Nratio,anchor,seed
    integer,intent(out) :: allign,debug,voltyp,domainw
    real(dp),intent(out) :: T,pres,epses,eps0,rsphere,spmyy,epsphere 
    real(dp),intent(out) :: sigma0,siges,Kw,cutoff,maxdr
    character(len=*), parameter :: paramsfile='gbcyl.in'
    integer, parameter :: input=20
    character(len=30),intent(out) :: file  
    integer,parameter :: parametercount=22
    integer :: readcount
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
    readcount=0
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
       read(unit=line,fmt=*,iostat=ios) string,file; readcount=readcount+1;
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
	 Nrelax=int(x+0.5); readcount=readcount+1;
       else if (string=='$Nprod') then
	 Nprod=int(x+0.5); readcount=readcount+1;
       else if (string=='$Nratio') then
	 Nratio=int(x+0.5); readcount=readcount+1;
       else if (string=='$T') then
	 T=x; readcount=readcount+1;
       else if (string=='$pres') then
	 pres=x; readcount=readcount+1;
       else if (string=='$anchor') then
	 anchor=int(x+0.5); readcount=readcount+1;
       else if (string=='$voltyp') then
         voltyp=int(x+0.5); readcount=readcount+1; 
       else if (string=='$seed') then
         seed=int(x+0.5); readcount=readcount+1;
       else if(string=='$epses')then 
         epses=x; readcount=readcount+1;
       else if(string=='$eps0')then
         eps0=x; readcount=readcount+1;
       else if(string=='$rsphere')then
         rsphere=x; readcount=readcount+1;
       else if(string=='$spmyy')then
         spmyy=x; readcount=readcount+1;
       else if(string=='$epsphere')then
         epsphere=x; readcount=readcount+1;
       else if(string=='$allign')then
         allign=int(x+0.5); readcount=readcount+1;
       else if(string=='$debug')then
         debug=int(x+0.5); readcount=readcount+1;
       else if(string=='$sigma0')then
         sigma0=x; readcount=readcount+1;
       else if(string=='$siges')then
         siges=x; readcount=readcount+1;
       else if(string=='$Kw')then
         Kw=x; readcount=readcount+1;
       else if(string=='$cutoff')then
         cutoff=x; readcount=readcount+1;
       else if(string=='$maxdr')then
         maxdr=x; readcount=readcount+1;
       else if(string=='$domainw')then
         domainw=int(x+0.5); readcount=readcount+1;
       else
         print '(A,A)','Unknown parameter',string
	 stop 'Parameter read in error'
       endif
     end if
       !!print '(A,A16,A,G13.6)','Read in parameter ',string,' value',x
    enddo
    if (readcount<parametercount) then 
      stop 'ReadParams: Some parameters are undefined'
    end if
    close(input)
  end subroutine ReadParams



subroutine writestateMPI(Rc,Lz,nall,GB,Xe,parray,np)
use particle, only : particledat
use nrtype
use mpi
implicit none
type(particledat), dimension(:), intent(in) :: parray
real(dp),intent(in) :: Rc,Lz
integer, intent(in) :: np,nall
integer,intent(in) :: GB,Xe
character(len=30) :: outfile='gbxecyl.out'
integer :: ios, outunit=19,i,astat,ntask,my_id,ierror
integer,dimension(:),allocatable :: help
integer :: s
!1. Kirjoitetaan partikkelimaarat ja sylinterin tiedot.
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ntask,ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierror)
  call MPI_REDUCE(np,s,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(my_id==0) then 
    write (*,*) s,nall
    if (s/=nall) stop 'nall/= what we are writing in to the file'
  end if   
  allocate(help(np),stat=astat)
  if(astat/=0) stop 'writestateMPI:allocation of array help failed'
  do i=1,np
    if(parray(i)%rod) then 
      help(i)=1
    else
      help(i)=0
    end if
  end do
  if(my_id==0) then 
    open(outunit,file=outfile,status='replace',position='append',&
         form='formatted',iostat=ios)
    if(ios/=0) then 
      write(*,*) 'writestateMPI: couldnt open file',outfile
      stop;
    end if
!    write(outunit,*) '$R:',Rc,'$Lz:',Lz
!    write(outunit,*) '$N:',nall,'$GB:',GB,'$Xe:',Xe
    partccount=nall
    GBcount=GB
    Xecount=Xe
    radius=Rc
    height=Lz
    write(outunit,nml=gbxecyl)
    write(outunit,*) '$x:'
  !2. Kirjoitetaan koordinaatit
    write(outunit,'(G12.6)',ADVANCE='YES') parray(1:np)%x
    close(outunit)
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  do i=1,ntask-1
    if (my_id==i) then
      open(outunit,file=outfile,status='old',position='append',action='write',&
           form='formatted',iostat=ios)
      if(ios/=0) then
        write(*,*) 'writestateMPI: rank',my_id,'couldnt open the file.'
        stop
      end if 
      write(outunit,'(G12.6)',ADVANCE='YES') parray(1:np)%x
      close(outunit)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end do
  do i=0,ntask-1
    if(my_id==i) then
      open(outunit,file=outfile,status='old',position='append',action='write',&
           form='formatted',iostat=ios)   
      if (my_id==0) write(outunit,*) '$y:'
      write(outunit,'(G12.6)',ADVANCE='YES') parray(1:np)%y
      close(outunit)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end do
  do i=0,ntask-1
    if(my_id==i) then
      open(outunit,file=outfile,status='old',position='append',action='write',&
           form='formatted',iostat=ios)   
      if (my_id==0) write(outunit,*) '$z:'
      write(outunit,'(G12.6)',ADVANCE='YES') parray(1:np)%z
      close(outunit)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end do
  do i=0,ntask-1
    if(my_id==i) then
      open(outunit,file=outfile,status='old',position='append',action='write',&
           form='formatted',iostat=ios)   
      if (my_id==0) write(outunit,*) '$rod:'
      write(outunit,'(G12.6)',ADVANCE='YES') help(1:np)
      close(outunit)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end do
  do i=0,ntask-1
    if(my_id==i) then
      open(outunit,file=outfile,status='old',position='append',action='write',&
           form='formatted',iostat=ios)   
      if (my_id==0) write(outunit,*) '$ux:'
      write(outunit,'(G12.6)',ADVANCE='YES') parray(1:np)%ux
      close(outunit)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end do
  do i=0,ntask-1
    if(my_id==i) then
      open(outunit,file=outfile,status='old',position='append',action='write',&
           form='formatted',iostat=ios)   
      if (my_id==0) write(outunit,*) '$uy:'
      write(outunit,'(G12.6)',ADVANCE='YES') parray(1:np)%uy
      close(outunit)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end do
  do i=0,ntask-1
    if(my_id==i) then
      open(outunit,file=outfile,status='old',position='append',action='write',&
           form='formatted',iostat=ios)   
      if (my_id==0) write(outunit,*) '$uz:'
      write(outunit,'(G12.6)',ADVANCE='YES') parray(1:np)%uz
      close(outunit)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end do
end subroutine writestateMPI

  
  
end module io
