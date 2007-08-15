
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
