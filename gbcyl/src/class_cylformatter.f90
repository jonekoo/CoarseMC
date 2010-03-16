module class_cylformatter
use nrtype
use particle
use utils
use class_parameter_writer
use class_parameterizer
use class_poly_box
implicit none
private

!! This module provides the following services:
!!
!! 1. readstate: Reading a geometry from a file. 
!! 2. findlast: Finding the position of last geometry in a file containing 
!!    multiple geometries.
!! 3. writestate: Writing a geometry to a file. 
!!
!! By geometry we mean in this module a cylindrical simulation box and 
!! ellipsoidal particles.  
!! cylformatter is meant to be used to hold the geometry file data.
!! 

public :: cylformatter
public :: readstate
public :: writestate
public :: init
!!public :: finalizeOutput
!!public :: io_write_parameters
public :: findlast
public :: beginmark
public :: endmark
public :: new_cylformatter
public :: delete

type cylformatter
  private
  character(len = 200) :: file = 'simdata.out'
  integer :: unit = 19
  character(len = 3) :: beginmark = '$R:'
  character(len = 3) :: endmark = 'EOF'
end type

interface init
  module procedure init_with, init_molfile
end interface

interface findlast
  module procedure general_findlast, cf_findlast
end interface

interface delete
  module procedure cf_delete
end interface

contains

function new_cylformatter(statefile) result(cf)
  character(len = *), intent(in) :: statefile
  type(cylformatter) :: cf
  integer :: ios
  cf%file = statefile
  cf%unit = 17
  open(file = cf%file, unit = cf%unit, status = 'old', action = 'readwrite', &
  iostat = ios)
  if (ios /= 0) then
    write(*, *) 'Could not open file:', cf%file
    stop
  end if
end function

subroutine cf_delete(cf)
  type(cylformatter), intent(inout) :: cf
  cf%file = ''
  close(cf%unit)
end subroutine

subroutine init_molfile(cf, statefile)
  type(cylformatter), intent(inout) :: cf
  character(len = *), intent(in) :: statefile
  cf%file = statefile
  call open_output(cf)
end subroutine

subroutine init_with(cf, reader)
  type(cylformatter), intent(inout) :: cf
  type(parameterizer), intent(in) :: reader
  call get_parameter(reader, 'molecule_file', cf%file)
  call open_output(cf)
end subroutine

subroutine io_write_parameters(cf, writer)
  type(cylformatter), intent(in) :: cf
  type(parameter_writer), intent(in) :: writer
  call write_parameter(writer, 'molecule_file', cf%file)
end subroutine

subroutine finalizeOutput(cf)
  type(cylformatter), intent(in) :: cf
  close(cf%unit) 
end subroutine finalizeOutput

subroutine open_output(cf)
  type(cylformatter), intent(in) :: cf
  integer :: ios
  open(UNIT = cf%unit, FILE = cf%file, & 
    & POSITION = 'APPEND', FORM = 'FORMATTED', IOSTAT = ios) 
  if (ios /= 0) then
    write (*, *) 'Could not open file ', cf%file, '.'
    stop;
  end if
end subroutine open_output

!! Writes the coordinates of @p particles and the dimensions of the 
!! cylindrical simulation cell to the file specified by @p cf.
!! 
!! @p cf the cylformatter object defining the file to write to.
!! @p particles the array of particles to write to the file.
!! @p n_particles number of particles in array @p particles.
!! @p radius the radius of the cylinder.
!! @p height the height of the cylinder.
!! 
subroutine writestate(cf, particles, n_particles, radius, height)
  type(cylformatter), intent(in) :: cf
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
  write(cf%unit, '(A3, 1X,' // fmt_char_dp() // ', 1X, A4, 1X, ' // &
  fmt_char_dp() // ')') '$R:', radius, '$Lz:', height
  write(cf%unit, '(A3, 1X, I7, 1X, A4, 1X, I7, 1X, A4, 1X, I7)') '$N:', &
  n_particles, '$GB:', GB, '$Xe:', Xe
  write(cf%unit, *) '$x:'
  write(cf%unit, '(' // fmt_char_dp() // ')') particles(1:n_particles)%x
  write(cf%unit, *) '$y:'
  write(cf%unit, '(' // fmt_char_dp() // ')') particles(1:n_particles)%y
  write(cf%unit, *) '$z:'
  write(cf%unit, '(' // fmt_char_dp() // ')') particles(1:n_particles)%z
  write(cf%unit, *) '$rod:'
  write(cf%unit, *) help(1:n_particles)
  write(cf%unit, *) '$ux:'
  write(cf%unit, '(' // fmt_char_dp() // ')') particles(1:n_particles)%ux
  write(cf%unit, *) '$uy:'
  write(cf%unit, '(' // fmt_char_dp() // ')') particles(1:n_particles)%uy
  write(cf%unit, *) '$uz:'
  write(cf%unit, '(' // fmt_char_dp() // ')') particles(1:n_particles)%uz
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

function beginmark(cf) result(mark)
  type(cylformatter), intent(in) :: cf
  character(len = 3) :: mark
  mark = cf%beginmark
end function

function endmark(cf) result(mark)
  type(cylformatter), intent(in) :: cf
  character(len = 3) :: mark
  mark = cf%endmark
end function endmark

subroutine readstate(cf, particles, n_particles, radius, height)
  type(cylformatter), intent(in) :: cf
  type(particledat), dimension(:), pointer :: particles
  integer, intent(out) :: n_particles
  real(dp), intent(out) :: radius, height
  integer :: astat
  integer :: i
  character(len = 3) :: charvar
  integer, dimension(:), allocatable :: help
  read(cf%unit, *) charvar, radius, charvar, height
  read(cf%unit, *) charvar, n_particles
  read(cf%unit, *) charvar
  allocate(particles(n_particles), help(n_particles), stat = astat)
  if (astat /= 0) then
    write(*,*) 'readstate: Virhe varattaessa muistia: particles, help'     
    stop
  end if   
  read(cf%unit, *) particles(1:n_particles)%x
  read(cf%unit, *) charvar
  read(cf%unit, *) particles(1:n_particles)%y
  read(cf%unit, *) charvar,particles(1:n_particles)%z
  read(cf%unit, *) charvar
  read(cf%unit, *) help(1:n_particles)
  read(cf%unit, *) charvar
  read(cf%unit, *) particles(1:n_particles)%ux
  read(cf%unit, *) charvar,particles(1:n_particles)%uy
  read(cf%unit, *) charvar,particles(1:n_particles)%uz
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
!! findlast and read_state when needed. This only aggregates those two and
!! file opening/closing.
!!
subroutine read_last_state(filename, begin, end, simbox, particles, n_particles)
  character(len = *), intent(in) :: filename
  character(len = *), intent(in) :: begin
  character(len = *), intent(in) :: end
  type(particledat), dimension(:), pointer :: particles
  integer, intent(out) :: n_particles
  integer, parameter :: readunit = 17   
  integer :: ios
  logical :: isfound
  type(poly_box), intent(out) :: simbox
  open(readunit, file = filename, status = 'old', iostat = ios)
  if (ios /= 0) then
    write(6, *) 'Error: Could not open file ', filename, '. Stopping.'
    stop 
  end if
  call findlast(readunit, begin, end, isfound)
  call read_state(readunit, simbox, particles, n_particles)
  close(readunit) 
end subroutine

subroutine cf_findlast(cf, isfound)
  type(cylformatter), intent(in) :: cf
  logical, intent(out), optional :: isfound
  call findlast(cf%unit, beginmark(cf), endmark(cf), isfound)
end subroutine


!! Finds last part of file that is limited by lines starting with characters 
!! @p begin and @p end. The position of readunit is moved to the beginning of
!! line where @p begin is found. 
!!
!! Rules and restrictions:
!! The @p readunit should be a real file and not e.g. standard input.
!!
!! A special case:
!! If @p end has the value 'EOF' end is not searched but end of file is used
!! instead.
!!
!! @p readunit the unit to which the file is connected to. 
!! @p begin the first characters in the beginning line.
!! @p end the first characters in the ending line.
!! @p isfound evaluates to .true. only if both begin and end are found in the 
!! file and begin is on a line preceding end. 
!!
subroutine general_findlast(readunit, begin, end, isfound)
  integer, intent(in) :: readunit
  character(len = *), intent(in) :: begin
  character(len = *), intent(in) :: end
  logical, intent(out) :: isfound
  character(len = max(len(begin), len(end))) :: line
  integer :: ios
  integer :: linenumber
  isfound = .false.
  ! size of line to max begin, end
  line = ''
  linenumber = 0
  !! Go to end of file
  do 
    linenumber = linenumber + 1
    read(readunit, iostat = ios , fmt = *) line
    if (ios < 0) then
      exit
    end if
  end do
  if('EOF' /= end) then
    write(*, *) 'Backspace until end of geometry found'
    call findbackwards(readunit, end, linenumber, isfound)
    if (.not. isfound) then
      return
    end if
  end if
  write(*, *) 'Backspace until beginning of geometry found'
  call findbackwards(readunit, begin, linenumber, isfound)
end subroutine
 
!! Moves the position specifier of file connected to @p readunit to the 
!! beginning of the line with first characters @p linetofind. Searching is
!! done backwards from the current position of file connected to @p readunit.
!! 
!! @p readunit unit connected to a file opened with reading capability
!! @p linetofind the first characters of the line to be found
!! @p linenumber - 1 is the maximum number of lines searched backwards.  
!! @p isfound evaluates to .true. if line is found and to .false. if not.
!!
subroutine findbackwards(readunit, linetofind, linenumber, isfound)
  integer, intent(in) :: readunit
  character(len = *), intent(in) :: linetofind
  integer, intent(inout) :: linenumber
  logical, intent(out) :: isfound
  character(len = len(linetofind)) :: line
  character(len = 30) :: lenchar
  isfound = .false.
  write(lenchar, *) len(line)
  do 
    backspace(readunit) 
    read(readunit, fmt = '(A' // trim(adjustl(lenchar)) // ')') line
    backspace(readunit)
    if (line == linetofind) then
      isfound = .true.
      return
    end if
    linenumber = linenumber - 1
    if (linenumber < 1) then
      return
    end if    
  end do
end subroutine

end module
