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
  public :: io_write_parameters
  public :: find_last
  public :: read_configuration

  character(len = 50), private, save :: molecule_data_file_ = 'simdata.out'
  integer, private, save :: molecule_data_unit_ = 19

  interface init
    module procedure init_old, init_with, init_molfile
  end interface

  contains

  subroutine init_molfile(statefile)
    character(len = *), intent(in) :: statefile
    molecule_data_file_ = statefile
    call open_output
  end subroutine

  subroutine init_old
    integer :: id
    character(len = 10) :: id_char
    integer :: ierr
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    write(id_char, '(I9)') id
    molecule_data_file_ = 'simdata.out' // '.' // trim(adjustl(id_char))  
    call open_output
  end subroutine

  subroutine init_with(reader)
    type(parameterizer), intent(in) :: reader
    call get_parameter(reader, 'molecule_file', molecule_data_file_)
    if ('' == molecule_data_file_) then
      call init
    else 
      call open_output
    end if 
  end subroutine

  subroutine io_write_parameters(writer)
    type(parameter_writer), intent(in) :: writer
    call write_parameter(writer, 'molecule_file', molecule_data_file_)
  end subroutine

  subroutine finalizeOutput()
    close(molecule_data_unit_) 
  end subroutine finalizeOutput

  subroutine open_output
    integer :: ios
    open(UNIT = molecule_data_unit_, FILE = molecule_data_file_, & 
      & POSITION = 'APPEND', FORM = 'FORMATTED', IOSTAT = ios) 
    if (ios /= 0) then
      write (*, *) 'Could not open file ', molecule_data_file_, '.'
      stop;
    end if
  end subroutine open_output

  ! Kirjoittaa tilan eli sylinterin s‰teen, korkeuden ja partikkelitaulukon
  ! tiedostoon  
  subroutine writestate(particles, n_particles, radius, height)
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
    write(molecule_data_unit_, '(' // fmt_char_dp() // ')') &
    particles(1:n_particles)%x
    write(molecule_data_unit_, *) '$y:'
    write(molecule_data_unit_, '(' // fmt_char_dp() // ')') &
    particles(1:n_particles)%y
    write(molecule_data_unit_, *) '$z:'
    write(molecule_data_unit_, '(' // fmt_char_dp() // ')') &
    particles(1:n_particles)%z
    write(molecule_data_unit_, *) '$rod:'
    write(molecule_data_unit_, *) help(1:n_particles)
    write(molecule_data_unit_, *) '$ux:'
    write(molecule_data_unit_, '(' // fmt_char_dp() // ')') &
    particles(1:n_particles)%ux
    write(molecule_data_unit_, *) '$uy:'
    write(molecule_data_unit_, '(' // fmt_char_dp() // ')') &
    particles(1:n_particles)%uy
    write(molecule_data_unit_, *) '$uz:'
    write(molecule_data_unit_, '(' // fmt_char_dp() // ')') &
    particles(1:n_particles)%uz
    deallocate(help)
  end subroutine writestate

  subroutine read_configuration(read_unit, simbox, particles, n_particles)
    integer, intent(in) :: read_unit
    type(poly_box), intent(out) :: simbox
    type(particledat), dimension(:), pointer :: particles
    integer, intent(out) :: n_particles
    real(dp) :: radius, height
    integer :: astat, i
    character(len = 3) :: charvar
    integer, dimension(:), allocatable :: help
    character(len = 200) :: box_string
    read(read_unit, *) charvar, radius, charvar, height
    read(read_unit, *) charvar, n_particles
    write(box_string, '(A, 3' // fmt_char_dp() // ', A)') 'cylinder', &
    2._dp * radius, 2._dp * radius, height, ' F F T'
    call create_box(simbox, box_string)
    read(read_unit, *) charvar
    allocate(particles(n_particles), help(n_particles), stat = astat)
    if (astat /= 0) then
      stop 'Error! Allocation of memory failed in read_configuration'
    end if   
    read(read_unit,*) particles(1:n_particles)%x
    read(read_unit,*) charvar
    read(read_unit,*) particles(1:n_particles)%y
    read(read_unit,*) charvar,particles(1:n_particles)%z
    read(read_unit,*) charvar
    read(read_unit,*) help(1:n_particles)
    read(read_unit,*) charvar
    read(read_unit,*) particles(1:n_particles)%ux
    read(read_unit,*) charvar,particles(1:n_particles)%uy
    read(read_unit,*) charvar,particles(1:n_particles)%uz
    close(read_unit) 
    do i = 1, n_particles
      if(help(i) == 1) then
        particles(i)%rod = .true.
      else
        particles(i)%rod = .false.
      end if
    end do
  end subroutine

  subroutine readstate(filename, particles, n_particles, radius, height)
    character(len = *), intent(in) :: filename
    type(particledat), dimension(:), pointer :: particles
    integer, intent(out) :: n_particles
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
    read(readunit,*) charvar, n_particles
    read(readunit,*) charvar
    !Varaa muistin taulukolle
    allocate(particles(n_particles), help(n_particles), stat = astat)
    if (astat /= 0) then
      write(*,*) 'readstate: Virhe varattaessa muistia: particles, help'
      stop
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
    do i = 1, n_particles
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
  subroutine read_state(readunit, simbox, particles, n_particles)
    integer, intent(in) :: readunit
    type(poly_box), intent(out) :: simbox
    type(particledat), dimension(:), pointer :: particles 
    integer, intent(out) :: n_particles
    character(len = 3) :: charvar
    character(len = 500) :: particle_string
    character(len = 200) :: box_string
    integer :: astat, i
    type(particledat) :: particle_read
    read(readunit, '(A200)') box_string
    call create_box(simbox, box_string)
    read(readunit, *) charvar, n_particles
    allocate(particles(n_particles), stat = astat)
    if (astat /= 0) then
      stop 'Error: Could not allocate memory for particle array. Stopping.'
    end if
    do i = 1, n_particles
      read(readunit, '(A500)') particle_string
      particle_read = create_particle(particle_string)
      particles(i) = particle_read 
    end do
  end subroutine

  !! Reads the last state of @p particles and @p simbox written in @p filename.
  !! 
  !! @p filename the name of the file to be read from
  !! @p simbox the simulation box to be read.
  !! @p particles the dynamic array to which particles will be read
  !! @p n_particles the number of particles in @p particles
  !!
  !! :NOTE: it may be better to discard this routine and use the two routines
  !! find_last and read_state when needed. This only aggregates those two and
  !! file opening/closing.
  !!
  subroutine read_last_state(filename, begin, end, simbox, particles, &
  n_particles)
    character(len = *), intent(in) :: filename
    character(len = *), intent(in) :: begin
    character(len = *), intent(in) :: end
    type(particledat), dimension(:), pointer :: particles
    integer, intent(out) :: n_particles
    integer, parameter :: readunit = 17   
    integer :: ios
    type(poly_box), intent(out) :: simbox
    open(readunit, file = filename, status = 'old', iostat = ios)
    if (ios /= 0) then
      write(6, *) 'Error: Could not open file ', filename, '. Stopping.'
      stop 
    end if
    call find_last(readunit, begin, end)
    call read_state(readunit, simbox, particles, n_particles)
    close(readunit) 
  end subroutine

  subroutine find_last(read_unit, begin, end)
    integer, intent(in) :: read_unit
    character(len = *), intent(in) :: begin
    character(len = *), intent(in) :: end
    character(len = 80) :: line
    integer :: ios
    ! size of line to max begin, end
    !! Go to end of file
    line = ''
    do 
      read(read_unit, iostat = ios ,fmt = '(A80)') line
      if(ios < 0) then
        exit
      end if
    end do
    !! Backspace until end of geometry found
    if('EOF' /= end) then
      do while(end /= trim(adjustl(line)))
        backspace(read_unit, iostat = ios)
        backspace(read_unit, iostat = ios)
        if (0 /= ios) then
          !! This doesn't make a difference if the file ends or there is some 
          !! other reading problem! 
          write(6, *) 'find_last: no end geometry line found in file ' 
          stop
        end if
        read(read_unit, *) line
      end do
    end if
    !! Backspace until beginning of geometry found
    do while(begin /= trim(adjustl(line)))
      backspace(read_unit, iostat = ios) !! Should it be used twice?
      backspace(read_unit, iostat = ios)
      if (0 /= ios) then 
        !! This doesn't make a difference if the file ends or there is some 
        !! other reading problem! 
        write(6, *) 'find_last: no begin geometry line found in file ' 
        stop
      end if
      read(read_unit, *) line
    end do
    backspace(read_unit, iostat = ios)
  end subroutine
 
end module
