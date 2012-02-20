module class_parameter_writer
use nrtype
use utils
implicit none
private

public :: new_parameter_writer
public :: delete
public :: parameter_writer
public :: writeparameter
public :: writecomment

type parameter_writer
  private
  integer :: unit
end type

interface writeparameter
  module procedure pw_writereal, pw_writestring, pw_writeint, pw_writelogical
end interface

interface delete
  module procedure pw_delete
end interface

interface writecomment
  module procedure pw_writecomment
end interface

contains

function new_parameter_writer(pwunit) result(pw)
  !character(len = *), intent(in) :: parameterfilename
  integer, intent(in) :: pwunit
  type(parameter_writer) :: pw
  pw%unit = pwunit
end function

subroutine pw_delete(p)
  type(parameter_writer), intent(inout) :: p
  !! If opened
  close(p%unit) 
end subroutine

subroutine pw_writecomment(p, comment)
  type(parameter_writer), intent(in) :: p
  character(len = *), intent(in) :: comment
  write(p%unit, *) comment
end subroutine

subroutine pw_writestring(p, key, value)
  type(parameter_writer), intent(in) :: p
  character(len = *), intent(in) :: key
  character(len = *), intent(in) :: value
  call pw_writeparameter(p, key, '"'//trim(adjustl(value))//'"')  
end subroutine

subroutine pw_writeparameter(p, key, value)
  type(parameter_writer), intent(in) :: p
  character(len = *), intent(in) :: key
  character(len = *), intent(in) :: value
  character(len = 30) :: keyfield
  character(len = 50) :: valuefield
  write(keyfield, *) "$" // adjustl(key)
  write(valuefield, '(A50)') value
  write(p%unit, '(A30, A50)') adjustl(keyfield), adjustl(valuefield)
end subroutine

subroutine pw_writereal(p, key, value)
  type(parameter_writer), intent(in) :: p
  character(len = *), intent(in) :: key
  real(dp), intent(in) :: value
  character(len = 50) :: valuestring
  write(valuestring, '(' // fmt_char_dp() // ')') value
  call pw_writeparameter(p, key, valuestring) 
end subroutine

subroutine pw_writeint(p, key, value)
  type(parameter_writer), intent(in) :: p 
  character(len = *), intent(in) :: key 
  integer, intent(in) :: value
  character(len = 50) :: valuestring
  write(valuestring, '(' // fmt_char_int() // ')') value 
  call pw_writeparameter(p, key, valuestring)
end subroutine

subroutine pw_writelogical(p, key, value)
  type(parameter_writer), intent(in) :: p
  character(len = *), intent(in) :: key
  logical, intent(in) :: value
  character(len = 50) :: valuestring
  write(valuestring, *) value
  call pw_writeparameter(p, key, valuestring)
end subroutine

end module

