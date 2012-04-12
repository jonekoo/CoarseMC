module class_parameterizer
use nrtype
use m_fileunit
implicit none
private

public :: parameterizer
public :: new_parameterizer
public :: getparameter
public :: delete

character(len=*), parameter :: logfile = 'parameterizer.log'

type parameterizer
  private
  integer :: unit
  character(len = 80), dimension(:, :), pointer :: parameters
  integer :: logunit = 6
end type

interface readstring
  module procedure parameterizer_readstring
end interface

interface new_parameterizer
  module procedure new_parameterizers !, new_parameterizeri
end interface

interface getparameter
  module procedure parameterizer_getstring, parameterizer_getinteger, parameterizer_getreal, parameterizer_getlogical
end interface

interface delete
  module procedure parameterizer_delete
end interface

contains

!> Returns a parameter reader trying to read parameters from @param 
!! parameterfile. Optionally @param logfile can be given in which case the
!! parameterizer will write information about parameters read into a file
!! instead of standard output. 
!! 
!! @param parameterfile the file parameters are read from.
!! @param logfile the optional logfile.
!!
function new_parameterizers(parameterfile, logfile, iostat) result(p)
  character(len = *), intent(in) :: parameterfile
  character(len = *), intent(in), optional :: logfile
  integer, optional, intent(out) :: iostat
  integer :: ios
  type(parameterizer) :: p
  integer :: logunit
  p%unit = fileunit_getfreeunit()
  open(file=parameterfile, unit=p%unit, status='OLD', action='READ', iostat=ios)
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

subroutine parameterizer_delete(p)
  type(parameterizer), intent(inout) :: p
  close(p%unit)
  if(p%logunit /= 6) close(p%logunit)
end subroutine

!! :TODO: Change searching from file to searching from table.
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

subroutine parameterizer_getstring(p, key, value)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  character(len = *), intent(inout) :: value
  character(len = 50) :: string
  call readstring(p, key, string)
  if ('' == string) then 
    write(p%logunit, *) 'Using default value ', key, ' = ', value
  else
    value = string
  end if
end subroutine

subroutine parameterizer_getinteger(p, key, value)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  integer, intent(inout) :: value
  character(len = 50) :: string
  integer :: ios
  call readstring(p, key, string)
  if ('' == string) then
    write(p%logunit, *) 'Using default value ', key, ' = ', value
  else
    read(string, fmt = *, iostat = ios) value
    if(0 /= ios) then
      call conversionwarning(p, key, string, 'integer')
    end if
  end if
end subroutine

subroutine parameterizer_getreal(p, key, value)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  real(dp), intent(inout) :: value
  character(len = 50) :: string
  integer :: ios
  call readstring(p, key, string)
  read(string, fmt = *, iostat = ios) value
  if ('' == string) then
    write(p%logunit, *) 'Using default value ', key, ' = ', value
  else
    if(0 /= ios) then
      call conversionwarning(p, key, string, 'real')
    end if
  end if
end subroutine

subroutine parameterizer_getlogical(p, key, value)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  logical, intent(inout) :: value
  character(len = 50) :: string
  integer :: ios
  call readstring(p, key, string)
  read(string, fmt = *, iostat = ios) value
  if ('' == string) then
    write(p%logunit, *) 'Using default value ', key, ' = ', value
  else
    if(0 /= ios) then
      call conversionwarning(p, key, string, 'logical')
    end if
  end if
end subroutine

subroutine conversionwarning(p, key, value, typestring)
  type(parameterizer), intent(in) :: p
  character(len = *), intent(in) :: key
  character(len = *), intent(in) :: value
  character(len = *), intent(in) :: typestring
  write(p%logunit, *) 'Warning: Could not convert ' // trim(adjustl(key)) // ' = ' // &
    trim(adjustl(value)) // ' to ' // trim(adjustl(typestring)) // '.'
end subroutine

end module
