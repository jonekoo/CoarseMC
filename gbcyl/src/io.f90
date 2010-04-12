module io
  use nrtype
  use particle
  use utils
  use mpi
  use class_parameter_writer
  use class_parameterizer
  use class_poly_box
  implicit none
  private

  public :: readstate
  public :: writestate
  public :: init
  public :: finalizeOutput
  public :: io_writeparameters
  public :: findlast
  public :: readconfiguration

  character(len = 50), private, save :: moleculedatafile = 'simdata.out'
  integer, private, save :: moleculedataunit = 19

  interface init
    module procedure initold, initwith, initmolfile
  end interface

  interface readstate
    module procedure readstatecyl, readstategen
  end interface

  contains

  subroutine initmolfile(statefile)
    character(len = *), intent(in) :: statefile
    moleculedatafile = statefile
    call openoutput
  end subroutine

  subroutine initold
    integer :: id
    character(len = 10) :: idchar
    integer :: ierr
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    write(idchar, '(I9)') id
    moleculedatafile = 'simdata.out' // '.' // trim(adjustl(idchar))  
    call openoutput
  end subroutine

  subroutine initwith(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'molecule_file', moleculedatafile)
    if ('' == moleculedatafile) then
      call init
    else 
      call openoutput
    end if 
  end subroutine

  subroutine io_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writeparameter(writer, 'molecule_file', moleculedatafile)
  end subroutine

  subroutine finalizeOutput()
    close(moleculedataunit) 
  end subroutine finalizeOutput

  subroutine openoutput
    integer :: ios
    open(UNIT = moleculedataunit, FILE = moleculedatafile, & 
      & POSITION = 'APPEND', FORM = 'FORMATTED', IOSTAT = ios) 
    if (ios /= 0) then
      write (*, *) 'Could not open file ', moleculedatafile, '.'
      stop;
    end if
  end subroutine openoutput

  ! Kirjoittaa tilan eli sylinterin s‰teen, korkeuden ja partikkelitaulukon
  ! tiedostoon  
  subroutine writestate(particles, nparticles, radius, height)
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: nparticles    
    real(dp), intent(in) :: radius
    real(dp), intent(in) :: height
    integer :: GB = 0, Xe = 0, astat, i
    integer, dimension(:), allocatable :: help
    allocate(help(nparticles), stat = astat)
    GB = 0
    Xe = 0
    do i = 1, nparticles
      if(particles(i)%rod) then 
        GB=GB+1
        help(i)=1
      else 
        Xe=Xe+1
        help(i)=0
      end if
    end do 
    write(moleculedataunit, *) '$R:',radius,'$Lz:',height
    write(moleculedataunit, *) '$N:',nparticles,'$GB:',GB,'$Xe:',Xe
    write(moleculedataunit, *) '$x:'
    write(moleculedataunit, '(' // fmt_char_dp() // ')') &
    particles(1:nparticles)%x
    write(moleculedataunit, *) '$y:'
    write(moleculedataunit, '(' // fmt_char_dp() // ')') &
    particles(1:nparticles)%y
    write(moleculedataunit, *) '$z:'
    write(moleculedataunit, '(' // fmt_char_dp() // ')') &
    particles(1:nparticles)%z
    write(moleculedataunit, *) '$rod:'
    write(moleculedataunit, *) help(1:nparticles)
    write(moleculedataunit, *) '$ux:'
    write(moleculedataunit, '(' // fmt_char_dp() // ')') &
    particles(1:nparticles)%ux
    write(moleculedataunit, *) '$uy:'
    write(moleculedataunit, '(' // fmt_char_dp() // ')') &
    particles(1:nparticles)%uy
    write(moleculedataunit, *) '$uz:'
    write(moleculedataunit, '(' // fmt_char_dp() // ')') &
    particles(1:nparticles)%uz
    deallocate(help)
  end subroutine writestate

  subroutine readconfiguration(readunit, simbox, particles, nparticles)
    integer, intent(in) :: readunit
    type(poly_box), intent(out) :: simbox
    type(particledat), dimension(:), pointer :: particles
    integer, intent(out) :: nparticles
    real(dp) :: radius, height
    integer :: astat, i
    character(len = 3) :: charvar
    integer, dimension(:), allocatable :: help
    character(len = 200) :: boxstring
    read(readunit, *) charvar, radius, charvar, height
    read(readunit, *) charvar, nparticles
    write(boxstring, '(A, 3' // fmt_char_dp() // ', A)') 'cylinder', &
    2._dp * radius, 2._dp * radius, height, ' F F T'
    call createbox(simbox, boxstring)
    read(readunit, *) charvar
    allocate(particles(nparticles), help(nparticles), stat = astat)
    if (astat /= 0) then
      stop 'Error! Allocation of memory failed in readconfiguration'
    end if   
    read(readunit,*) particles(1:nparticles)%x
    read(readunit,*) charvar
    read(readunit,*) particles(1:nparticles)%y
    read(readunit,*) charvar,particles(1:nparticles)%z
    read(readunit,*) charvar
    read(readunit,*) help(1:nparticles)
    read(readunit,*) charvar
    read(readunit,*) particles(1:nparticles)%ux
    read(readunit,*) charvar,particles(1:nparticles)%uy
    read(readunit,*) charvar,particles(1:nparticles)%uz
    close(readunit) 
    do i = 1, nparticles
      if(help(i) == 1) then
        particles(i)%rod = .true.
      else
        particles(i)%rod = .false.
      end if
    end do
  end subroutine

  subroutine readstatecyl(filename, particles, nparticles, radius, height)
    character(len = *), intent(in) :: filename
    type(particledat), dimension(:), pointer :: particles
    integer, intent(out) :: nparticles
    real(dp), intent(out) :: radius, height
    integer,parameter :: readunit = 17   
    integer :: ios, astat, i
    character(len = 3) :: charvar
    integer, dimension(:), allocatable :: help
    !Koittaa avata annetun tiedoston
    open(readunit, file = filename, status = 'old', iostat = ios)
    if (ios /= 0) then
      write(*,*) 'Tiedoston ', filename, &
      ' avaaminen ep‰onnistui, virhekoodi: ', ios
      write(*,*) 'Ohjelman suoritus keskeytyy.'
      stop;
    end if
    !Luetaan hiukkasten lukum‰‰r‰ tiedostosta
    read(readunit,*) charvar, radius, charvar, height
    read(readunit,*) charvar, nparticles
    read(readunit,*) charvar
    !Varaa muistin taulukolle
    allocate(particles(nparticles), help(nparticles), stat = astat)
    if (astat /= 0) then
      write(*,*) 'readstate: Virhe varattaessa muistia: particles, help'
      stop
    end if   
    !Luetaan  hiukkasten x-koordinaatit
    read(readunit,*) particles(1:nparticles)%x
    read(readunit,*) charvar
    !Luetaan y-koordinaatit
    read(readunit,*) particles(1:nparticles)%y
    read(readunit,*) charvar,particles(1:nparticles)%z
    !Luetaan rod-tieto
    read(readunit,*) charvar
    read(readunit,*) help(1:nparticles)
    !Luetaan orientaatiovektoreiden komponentit
    read(readunit,*) charvar
    read(readunit,*) particles(1:nparticles)%ux
    read(readunit,*) charvar,particles(1:nparticles)%uy
    read(readunit,*) charvar,particles(1:nparticles)%uz
    close(readunit) 
    do i = 1, nparticles
      if(help(i) == 1) then
        particles(i)%rod = .true.
      else
        particles(i)%rod = .false.
      end if
    end do
  end subroutine

  !! Move this to another module. Think about how the serialization/writing of
  !! states should be done. The objects should be responsible of choosing the
  !! necessary variables to write on disk but the format should be independent
  !! of the object itself. So the mechanism could be very similar to the 
  !! parameterizer/parameter_reader mechanism. We could however impose an 
  !! additional rule that the data is written and read sequentially and the 
  !! order is dictated by the object to be written.
  !!
  subroutine readstategen(readunit, simbox, particles, nparticles)
    integer, intent(in) :: readunit
    type(poly_box), intent(out) :: simbox
    type(particledat), dimension(:), pointer :: particles 
    integer, intent(out) :: nparticles
    character(len = 3) :: charvar
    character(len = 500) :: particlestring
    character(len = 200) :: boxstring
    integer :: astat, i
    type(particledat) :: particleread
    read(readunit, '(A200)') boxstring
    call createbox(simbox, boxstring)
    read(readunit, *) charvar, nparticles
    allocate(particles(nparticles), stat = astat)
    if (astat /= 0) then
      stop 'Error: Could not allocate memory for particle array. Stopping.'
    end if
    do i = 1, nparticles
      read(readunit, '(A500)') particlestring
      particleread = createparticle(particlestring)
      particles(i) = particleread 
    end do
  end subroutine

  !! Reads the last state of @p particles and @p simbox written in @p filename.
  !! 
  !! @p filename the name of the file to be read from
  !! @p simbox the simulation box to be read.
  !! @p particles the dynamic array to which particles will be read
  !! @p nparticles the number of particles in @p particles
  !!
  !! :NOTE: it may be better to discard this routine and use the two routines
  !! findlast and readstate when needed. This only aggregates those two and
  !! file opening/closing.
  !!
  subroutine readlaststate(filename, begin, end, simbox, particles, &
  nparticles)
    character(len = *), intent(in) :: filename
    character(len = *), intent(in) :: begin
    character(len = *), intent(in) :: end
    type(particledat), dimension(:), pointer :: particles
    integer, intent(out) :: nparticles
    integer, parameter :: readunit = 17   
    integer :: ios
    type(poly_box), intent(out) :: simbox
    open(readunit, file = filename, status = 'old', iostat = ios)
    if (ios /= 0) then
      write(6, *) 'Error: Could not open file ', filename, '. Stopping.'
      stop 
    end if
    call findlast(readunit, begin, end)
    call readstate(readunit, simbox, particles, nparticles)
    close(readunit) 
  end subroutine

  subroutine findlast(readunit, begin, end)
    integer, intent(in) :: readunit
    character(len = *), intent(in) :: begin
    character(len = *), intent(in) :: end
    character(len = 80) :: line
    integer :: ios
    ! size of line to max begin, end
    !! Go to end of file
    line = ''
    do 
      read(readunit, iostat = ios ,fmt = '(A80)') line
      if(ios < 0) then
        exit
      end if
    end do
    !! Backspace until end of geometry found
    if('EOF' /= end) then
      do while(end /= trim(adjustl(line)))
        backspace(readunit, iostat = ios)
        backspace(readunit, iostat = ios)
        if (0 /= ios) then
          !! This doesn't make a difference if the file ends or there is some 
          !! other reading problem! 
          write(6, *) 'findlast: no end geometry line found in file ' 
          stop
        end if
        read(readunit, *) line
      end do
    end if
    !! Backspace until beginning of geometry found
    do while(begin /= trim(adjustl(line)))
      backspace(readunit, iostat = ios) !! Should it be used twice?
      backspace(readunit, iostat = ios)
      if (0 /= ios) then 
        !! This doesn't make a difference if the file ends or there is some 
        !! other reading problem! 
        write(6, *) 'findlast: no begin geometry line found in file ' 
        stop
      end if
      read(readunit, *) line
    end do
    backspace(readunit, iostat = ios)
  end subroutine
 
end module
