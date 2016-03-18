program readjsonlines
  use json_module
  use iso_fortran_env, only: error_unit
  implicit none
  character(len=:), allocatable :: line
  integer, parameter :: u_json = 11
  type(json_file) :: json_f
  integer :: i
  open(unit=u_json, file='test.json', action='READ')

  i = 0
  do while(ReadLine(u_json, line))
     !write(*, '(A)') line
     i = i + 1
     write(*, *) 'loading JSON from line ', i
     call json_f%load_from_string(line)
     call json_f%destroy()
  end do

  close(unit=u_json)
contains

  ! http://stackoverflow.com/a/21953596
  function ReadLine(aunit, InLine, trimmed) result(OK)
    integer, intent(IN) :: aunit
    character(LEN=:), allocatable, optional, intent(out) :: InLine
    logical, intent(in), optional :: trimmed
    integer, parameter :: line_buf_len= 1024*4
    character(LEN=line_buf_len) :: InS
    logical :: OK, set
    integer status, size
    
    OK = .false.
    set = .true.
    do
       read (aunit,'(a)',advance='NO',iostat=status, size=size) InS
       OK = .not. IS_IOSTAT_END(status)
       if (.not. OK) return
       if (IS_IOSTAT_END(status)) return
       if (present(InLine)) then
          if (set) then
             InLine = InS(1:size)
             set=.false.
          else
             InLine = InLine // InS(1:size)
          end if
       end if
       if (IS_IOSTAT_EOR(status)) exit
    end do
    if (present(trimmed) .and. present(InLine)) then
       if (trimmed) InLine = trim(adjustl(InLine))
    end if
    
  end function ReadLine
end program readjsonlines
