!> Implements functions for reading parameters from a file. Conforms to
!! the same format as the module class_parameter_writer. 
module class_parameterizer
use num_kind
use m_fileunit
implicit none
private

public :: parameterizer
public :: new_parameterizer
public :: getparameter
public :: delete

character(len=*), parameter :: logfile = 'parameterizer.log'

!> Holds the input unit for reading parameters and the output unit for
!! logging what was read. 
type parameterizer
  private
  integer :: unit
  character(len = 80), dimension(:, :), pointer :: parameters
  integer :: logunit = 6
end type

!> Generic interface for reading a parameter value as a string from the
!! input, so that minimal changes are needed to change input source to
!! e.g. a list of strings.
interface readstring
  module procedure parameterizer_readstring
end interface

interface new_parameterizer
  module procedure new_parameterizers
end interface

!> Interface for reading parameters.
interface getparameter
  module procedure parameterizer_getstring, parameterizer_getinteger, &
       parameterizer_getreal, parameterizer_getlogical
end interface

interface delete
  module procedure parameterizer_delete
end interface

contains

!> Returns a parameterizer which reads parameters from @p parameterfile.
!! Optionally @p logfile can be given in which case the parameterizer
!! will write information about parameters read into a file instead of
!! standard output. An error has occurred, if @p iostat /= 0 at return.
function new_parameterizers(parameterfile, logfile, iostat) result(p)
  character(len = *), intent(in) :: parameterfile
  character(len = *), intent(in), optional :: logfile
  integer, optional, intent(out) :: iostat
  integer :: ios
  type(parameterizer) :: p
  integer :: logunit
  p%unit = fileunit_getfreeunit()
  open(file=parameterfile, unit=p%unit, status='OLD', action='READ', & 
       iostat=ios)
  if (present(iostat)) iostat=ios
  if (0/=ios) then 
    write(*, *) 'new_parameterizers: Failed opening ', parameterfile
  else
    rewind p%unit
    if(present(logfile)) then
      logunit = fileunit_getfreeunit()
      open(file=logfile, unit=logunit, position='APPEND', status='UNKNOWN', &
      action='WRITE', iostat=ios)
      if (ios == 0) then
        p%logunit = logunit
      else
        write(p%logunit, *) 'Warning: Opening logfile ', logfile, ' failed.'
      end if
    end if
  end if
end function

!> Closes the output and input units of the parameterizer @p p.
subroutine parameterizer_delete(p)
  type(parameterizer), intent(inout) :: p
  close(p%unit)
  if(p%logunit /= 6) close(p%logunit)
end subroutine

!> Searches the input unit of @p p for a parameter with @p key and
!! @p value. If @p key is not found @p value is an empty string and
!! a warning is written to the logunit of @p p.
!!
!! @todo@ Change searching from file to searching from table.
!! @todo@ Return error code when key is not found.
!!
subroutine parameterizer_readstring(p, key, value)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  character(len = *), intent(inout) :: value
  character(len = 80) :: line
  character(len = 30) :: field1
  character(len = 50) :: field2 
  integer :: i
  integer :: ios
  value = ''
  rewind(p%unit, iostat = ios) 
  if(0 /= ios) then
    stop 'Can not rewind parameter file.'
  else
    i = 0
    do
       i = i + 1
       ! First read in line in whole.
       read(p%unit, fmt = '(A80)', iostat = ios) line
       ! Check for probable end-of-file
       if (ios < 0) then
         write(p%logunit, *) 'Warning! Parameter ', trim(key), ' not found!'
         exit
       end if
       ! If not parameter line, cycle
       if (line(1:1) /= '$') cycle
       ! Attempt to parse line into command and value
       ! using internal formatted read
       read(unit = line, fmt=*, iostat = ios) field1, field2
       if (ios > 0) then
          write(*, *) 'ERROR: Data read in error on line ', i
          stop 'Invalid parameter file'
       endif
       ! From here on line should be in variable format
       ! Look for variable-defining strings
       if(field1 == '$' // trim(key)) then
         value = trim(adjustl(field2))
         write(p%logunit, *) 'Read in parameter ', trim(key), ' = ', value
         exit
       end if
     end do
  end if
end subroutine


!> Reads a string parameter with @p key to @p value from the input
!! specified by @p p. @p found is .true. if the parameter was found.
subroutine parameterizer_getstring(p, key, value, found)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  character(len = *), intent(inout) :: value
  logical, optional, intent(out) :: found
  character(len = 50) :: string
  call readstring(p, key, string)
  if ('' == string) then 
    write(p%logunit, *) 'Using default value ', key, ' = ', value
    if (present(found)) found = .false.
  else
    value = string
    if (present(found)) found = .true.
  end if
end subroutine

!> Reads an integer parameter with @p key to @p value from the input
!! specified by @p p. @p found is .true. if the parameter was found.
subroutine parameterizer_getinteger(p, key, value, found)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  integer, intent(inout) :: value
  logical, optional, intent(out) :: found
  character(len = 50) :: string
  integer :: ios
  call readstring(p, key, string)
  read(string, fmt = *, iostat = ios) value
  if ('' == string) then
    write(p%logunit, *) 'Using default value ', key, ' = ', value
    if (present(found)) found = .false.
  else
    if(0 /= ios) then
      call conversionwarning(p, key, string, 'integer')
    end if
    if (present(found)) found = .true.
  end if
end subroutine

!> Reads a real type parameter with @p key to @p value from the input
!! specified by @p p. @p found is .true. if the parameter was found.
subroutine parameterizer_getreal(p, key, value, found)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  real(dp), intent(inout) :: value
  logical, optional, intent(out) :: found
  character(len = 50) :: string
  integer :: ios
  call readstring(p, key, string)
  read(string, fmt = *, iostat = ios) value
  if ('' == string) then
    write(p%logunit, *) 'Using default value ', key, ' = ', value
    if (present(found)) found = .false.
  else
    if(0 /= ios) then
      call conversionwarning(p, key, string, 'real')
    end if
    if (present(found)) found = .true.
  end if
end subroutine

!> Reads a logical-type parameter with @p key to @p value from the input
!! specified by @p p. @p found is .true. if the parameter was found.
subroutine parameterizer_getlogical(p, key, value, found)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  logical, intent(inout) :: value
  logical, optional, intent(out) :: found
  character(len = 50) :: string
  integer :: ios
  call readstring(p, key, string)
  read(string, fmt = *, iostat = ios) value
  if ('' == string) then
    write(p%logunit, *) 'Using default value ', key, ' = ', value
    if (present(found)) found = .false.
  else
    if(0 /= ios) then
      call conversionwarning(p, key, string, 'logical')
    end if
    if (present(found)) found = .true.
  end if
end subroutine

!> Writes a warning to the logunit of @p p if @p value of the parameter
!! with @p key could not be converted to the format given in
!! @p typestring.
subroutine conversionwarning(p, key, value, typestring)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  character(len = *), intent(in) :: value
  character(len = *), intent(in) :: typestring
  write(p%logunit, *) 'Warning: Could not convert ' // trim(adjustl(key)) // ' = ' // &
    trim(adjustl(value)) // ' to ' // trim(adjustl(typestring)) // '.'
end subroutine

end module
