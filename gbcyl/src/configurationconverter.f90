program configurationconverter
use class_poly_box
use particle
use class_factory
use m_fileunit
implicit none

type cylformatter
  character(len = 200) :: file = ''
  integer :: inunit = 5
  integer :: outunit = 6
  character(len = 3) :: beginmark = '$R:'
  character(len = 3) :: endmark = 'EOF'
end type


type(poly_box) :: simbox
type(particledat), dimension(:), pointer :: particles
type(cylformatter) :: reader
type(factory) :: writer
integer :: ios

reader = new_cylformatter()
!! Read configuration in with cylformatter
do 
  call cf_readstate(reader, simbox, particles, ios)
  if (ios /= 0) exit
  !! Print configuration out with factory
  !writer = new_factory()
  call factory_writestate(writer, 6, simbox, particles)
end do

contains

function new_cylformatter(statefile) result(cf)
  character(len = *), intent(in), optional :: statefile
  type(cylformatter) :: cf
  integer :: ios
  if (present(statefile)) then
    cf%file = statefile
    cf%inunit = fileunit_getfreeunit()
    cf%outunit= cf%inunit
    open(file = cf%file, unit = cf%inunit, action = 'READWRITE', &
    iostat = ios)
    if (ios /= 0) then
      write(*, *) 'Could not open file: ', cf%file
      stop
    end if
  end if
end function

!! Reads the cylinder dimensions and particle coordinates from a unit pointed 
!! to by the cylformatter object.
!! 
!! @p cf the cylformatter object that holds the unit to be read from.
!! @p simbox at return this is assigned to the cylinder object read from 
!!    the input unit.
!! @p particles the particles read from the input unit
!! @p iostatus the status of the read operation at return. Negative value
!! means an end of file condition and a positive value indicates an error as
!! dictated in the Fortran 90/95 standard.
!!
subroutine cf_readstate(cf, simbox, particles, iostatus)
  type(cylformatter), intent(in) :: cf
  type(poly_box), intent(out) :: simbox
  type(particledat), dimension(:), pointer :: particles
  integer, intent(out) :: iostatus
  integer :: nparticles
  real(dp) :: radius, height
  integer :: astat
  integer :: i
  character(len = 3) :: charvar
  character(len = 4) :: charvar2
  integer, dimension(:), allocatable :: help
  read(cf%inunit, fmt=*, iostat = iostatus) charvar, radius, charvar2, height!, anothercharvar
  if (iostatus < 0) return
  simbox = new_cylinder(2._dp * radius, height)
  read(cf%inunit, *, iostat = iostatus) charvar, nparticles
  read(cf%inunit, *, iostat = iostatus) charvar
  if(associated(particles)) then
    if (size(particles) /= nparticles) then
      deallocate(particles)
      particles=>NULL()
    end if
  end if
  if( .not. associated(particles)) then
    allocate(particles(nparticles), stat = astat)
    if (astat /= 0) then
      stop 'cf_readstate: Memory allocation for particles failed.'
    end if
  end if   
  allocate(help(nparticles), stat = astat)
  if (astat /= 0) then
    stop 'cf_readstate: Memory allocation for help failed.'
  end if
  read(cf%inunit, *, iostat = iostatus) particles(1:nparticles)%x
  read(cf%inunit, *, iostat = iostatus) charvar
  read(cf%inunit, *, iostat = iostatus) particles(1:nparticles)%y
  read(cf%inunit, *, iostat = iostatus) charvar,particles(1:nparticles)%z
  read(cf%inunit, *, iostat = iostatus) charvar
  read(cf%inunit, *, iostat = iostatus) help(1:nparticles)
  read(cf%inunit, *, iostat = iostatus) charvar
  read(cf%inunit, *, iostat = iostatus) particles(1:nparticles)%ux
  read(cf%inunit, *, iostat = iostatus) charvar,particles(1:nparticles)%uy
  read(cf%inunit, *, iostat = iostatus) charvar,particles(1:nparticles)%uz
  do i = 1, nparticles
    if(help(i) == 1) then
      particles(i)%rod = .true.
    else
      particles(i)%rod = .false.
    end if
  end do
end subroutine

end program
