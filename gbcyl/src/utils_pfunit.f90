module utils_pfunit
  use pfunit
  use utils
  use m_fileunit, only: fileunit_getfreeunit
  implicit none


contains

  subroutine test_readstr
    character(len=:), allocatable :: str
    character(len=12), parameter :: ref = 'first' // new_line(ref) // 'second'
    integer :: unit
    integer :: ios
    character(len=80) :: filename
    unit = fileunit_getfreeunit()
    open(unit=unit, file='test_readstr.txt', status='replace', &
         action='readwrite', access='stream', form='formatted')
    write(unit, '(A)', advance='no') ref
    rewind(unit)
    call readstr(unit, str, ios)
    call assertEqual(ref // new_line(ref), str, &
         'Written and read strings do not match.')
    close(unit) !, status='delete')
  end subroutine test_readstr

  subroutine test_readlines
    character(len=:), allocatable :: str
    character(len=13), parameter :: ref = 'qwerty'
    integer :: unit
    integer :: ios
    character(len=80) :: filename
    type(str_wrapper), allocatable :: lines(:)
    character(len=:), allocatable :: temp
    unit = fileunit_getfreeunit()
    open(unit=unit, status='scratch', action='readwrite')
    write(unit, '(A)') 'first'
    write(unit, '(A)') 'second'
    rewind(unit)
    call readstr(unit, temp, ios)
    rewind(unit)
    
    call readlines(unit, lines, ios)
    call assertEqual(2, size(lines), 'Number of lines does not match')
    if (size(lines) > 0) then
       call assertEqual('first', lines(1)%c, 'First line read does not match')
    end if
    if (size(lines) > 1) then
       call assertEqual('second', lines(2)%c, 'Second line read does not match')
    end if
    close(unit)
  end subroutine test_readlines
  
end module utils_pfunit
