program test_json_leak
  use json_module
  implicit none

  type(json_value), pointer :: json_val
  type(json_file) :: json
  character(kind=CK, len=:), allocatable :: str
  integer :: i, u_output
  
  call json_initialize()
  open(unit=u_output, file='json_leak.json')
  do i = 1, 2
     call json%load_file(filename='input-0.json')
     call json%get('', json_val)
     call json_print_to_string(json_val, str)
     write(u_output, advance='no', fmt='(A)') str
     call json%destroy()
  end do
  close(u_output)
end program test_json_leak
