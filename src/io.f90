module io
  use num_kind, only: dp
  use particle, only: particledat
  implicit none


  character(len = 15), private, save :: molecule_data_file_
  integer, private, save :: molecule_data_unit_
  logical, private, save :: debug_ = .false. 
  namelist /io_nml/ molecule_data_file_, molecule_data_unit_, debug_
  private :: io_nml



  contains

  

  subroutine init
    implicit none
    molecule_data_file_ = "simdata.out"
    molecule_data_unit_ = 19
    call open_output
  end subroutine init



  subroutine finalizeOutput()
    implicit none
    close(molecule_data_unit_) 
  end subroutine finalizeOutput



  !! Saves the module state to @p write_unit
  !!
  !! Precondition: Module state must be initialized.
  !! 
  subroutine save_state(write_unit)
    implicit none
    integer, intent(in) :: write_unit
    write(write_unit, NML = io_nml)
  end subroutine save_state



  !! Restores the module state from @p read_unit
  !!
  !! Precondition: module state must not be initialized by other means
  !! Postcondition: module state is initialized. 
  !!  
  subroutine load_state(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    read(read_unit, NML = io_nml)
    call open_output
  end subroutine load_state



  subroutine open_output
    implicit none
    integer :: ios
    open(UNIT = molecule_data_unit_, FILE = molecule_data_file_, & 
      & POSITION = 'APPEND', FORM = 'FORMATTED', IOSTAT = ios) 
    if (ios /= 0) then
      write (*, *) 'Could not open file ', molecule_data_file_, '.'
      stop;
    end if
  end subroutine open_output



  !Kirjoittaa tilan eli sylinterin s‰teen, korkeuden ja partikkelitaulukon
  !tiedostoon  
  subroutine writestate(particles, n_particles, radius, height)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles    
    real(dp), intent(in) :: radius
    real(dp), intent(in) :: height
    integer :: GB = 0, Xe = 0, astat, i
    integer, dimension(:), allocatable :: help
    allocate(help(n_particles), stat = astat)
    GB = 0
    Xe = 0
    do i = 1, n_particles
      if(particles(i)%rod) then 
        GB=GB+1
        help(i)=1
      else 
        Xe=Xe+1
        help(i)=0
      end if
    end do 
    write(molecule_data_unit_, *) '$R:',radius,'$Lz:',height
    write(molecule_data_unit_, *) '$N:',n_particles,'$GB:',GB,'$Xe:',Xe
    write(molecule_data_unit_, *) '$x:'
    write(molecule_data_unit_, '(G12.6)') particles(1:n_particles)%x
    write(molecule_data_unit_, *) '$y:'
    write(molecule_data_unit_, '(G12.6)') particles(1:n_particles)%y
    write(molecule_data_unit_, *) '$z:'
    write(molecule_data_unit_, '(G12.6)') particles(1:n_particles)%z
    write(molecule_data_unit_, *) '$rod:'
    write(molecule_data_unit_, *) help(1:n_particles)
    write(molecule_data_unit_, *) '$ux:'
    write(molecule_data_unit_, '(G12.6)') particles(1:n_particles)%ux
    write(molecule_data_unit_, *) '$uy:'
    write(molecule_data_unit_, '(G12.6)') particles(1:n_particles)%uy
    write(molecule_data_unit_, *) '$uz:'
    write(molecule_data_unit_, '(G12.6)') particles(1:n_particles)%uz
    deallocate(help)
  end subroutine writestate



  subroutine readstate(filename, particles, n_particles, radius, height)
    implicit none
    character(len=*), intent(in) :: filename
    type(particledat), dimension(:), pointer :: particles
    integer, intent(out) :: n_particles
    real(dp), intent(out) :: radius,height
    integer,parameter :: readunit=17   
    integer :: ios,astat,i
    character(len=3) :: charvar
    integer,dimension(:),allocatable :: help
    radius=0.
    height=0.  
    !Koittaa avata annetun tiedoston
    open(readunit,file=filename,status='old',iostat=ios)
    if (ios/=0) then
      write(*,*) 'Tiedoston ',filename,' avaaminen ep‰onnistui.'
      write(*,*) 'Ohjelman suoritus keskeytyy.'
      stop;
    end if
    !Luetaan hiukkasten lukum‰‰r‰ tiedostosta
    read(readunit,*) charvar,radius,charvar,height
    read(readunit,*) charvar,n_particles
    read(readunit,*) charvar
    !Varaa muistin taulukolle
    allocate(particles(n_particles),help(n_particles),stat=astat)
    if (astat/=0) then
      write(*,*) 'readstate:Virhe varattaessa muistia: particles,help'
      stop;
    end if   
    !Luetaan  hiukkasten x-koordinaatit
    read(readunit,*) particles(1:n_particles)%x
    read(readunit,*) charvar
    !Luetaan y-koordinaatit
    read(readunit,*) particles(1:n_particles)%y
    read(readunit,*) charvar,particles(1:n_particles)%z
    !Luetaan rod-tieto
    read(readunit,*) charvar
    read(readunit,*) help(1:n_particles)
    !Luetaan orientaatiovektoreiden komponentit
    read(readunit,*) charvar
    read(readunit,*) particles(1:n_particles)%ux
    read(readunit,*) charvar,particles(1:n_particles)%uy
    read(readunit,*) charvar,particles(1:n_particles)%uz
    close(readunit) 
    do i=1,n_particles
      if(help(i)==1) then
        particles(i)%rod=.true.
      else
        particles(i)%rod=.false.
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



  !Aliohjelma parameterien lukemiseen tiedostosta paramsfile.
  !Tekij‰: Juho Lintuvuori
  !muutokset: Jouni Karjalainen
  subroutine ReadParams(file, Nrelax, Nprod, Nratio, T, pres, anchor, voltyp,&
                        Kw, seed, epses, eps0, rsphere, spmyy, epsphere, &
                        sigma0, siges, allign, debug, adjusttype, moveratio,&
                        scalingratio, maxtranslation, maxrotation)
    implicit none
    ! Alkuper√§inen aliohjelma Antti Kurosen k√§sialaa; kurssilta johdatus
    ! atomistisiin simulaatioihin

    integer,intent(out) :: Nrelax,Nprod,Nratio,anchor,seed, adjusttype
    integer,intent(out) :: allign,debug,voltyp
    real(dp),intent(out) :: T,pres,epses,eps0,rsphere,spmyy,epsphere 
    real(dp),intent(out) :: sigma0,siges,Kw
    real(dp), intent(out) :: moveratio, scalingratio 
    real(dp), intent(out) :: maxtranslation, maxrotation
    character(len=*), parameter :: paramsfile='gbcyl.in'
    integer, parameter :: input=20
    character(len = 78), intent(out) :: file  

    integer :: i, ios
    character(len=100) :: line
    character(len=80) :: string
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
       else if(string=='$allign')then
         allign=int(x+0.5);
       else if(string=='$debug')then
         debug=int(x+0.5);
       else if(string=='$sigma0')then
         sigma0=x;
       else if(string=='$siges')then
         siges=x;
       else if(string=='$Kw')then
         Kw=x;
       else if(string=='$adjusttype')then
         adjusttype=int(x+0.5);
       else if(string=='$moveratio')then
         moveratio=x;
       else if(string=='$scalingratio')then
         scalingratio=x;
       else if(string=='$maxtranslation')then
         maxtranslation=x;
       else if(string=='$maxrotation') then
         maxrotation=x;
       else
         print '(A,A)','Unknown parameter',string
         stop 'Parameter read in error'
       endif
     end if
     if(debug_) print '(A,A16,A,G13.6)','Read in parameter ',string,' value',x
    enddo
    close(input)
  end subroutine ReadParams

  

end module io
