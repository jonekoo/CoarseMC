program example_jsonfortran
  use iso_fortran_env
  use json_module

  type(json_file) :: json
  logical :: found
  !integer :: i,j,k
  character(kind=CK, len=:), allocatable :: interaction_type
  type(json_value), pointer :: p, child
  integer(INT32) :: i
  character(kind=CK, len=1) :: ic
  ! initialize the module
  call json_initialize()
        
  ! read the file
  call json%load_file(filename = 'test1.json')

  ! print the file to the console
  !call json%print_file()
  
  ! extract data from the file
  ! [found can be used to check if the data was really there]
  
  ! Three ways to list the types of interactions present:
  
  call json%get('interactions', p, found)
  
  write(output_unit, '(/, A)') '1) With callback'
  call json_get(p, print_type)
  
  write(output_unit, '(/, A)') '2) The ugly but simple way.'

  do i = 1, json_count(p)
     write(ic, '(I1)') i
     call json%get('interactions('// ic //').type', interaction_type, found)
     if ( .not. found ) then
        stop 1
     else
        write(output_unit, *) interaction_type
     end if
  end do
  
  write(output_unit, '(/, A)') '3) Iterating through the children manually.'
  do i = 1, json_count(p) 
     call json_get_child(p, i, child)
     call json_get(child, 'type', interaction_type, found)
     if (.not. found) then
        stop 1
     end if
     write(output_unit, *) interaction_type
  end do
  
  ! clean up
  call json%destroy()
  if (json_failed()) stop 1

contains
  subroutine print_type(element, i, count)
    type(json_value), pointer, intent(in) :: element
    integer(INT32), intent(in) :: i, count
    call json_get(element, 'type', interaction_type, found)
    write(output_unit, *) interaction_type
  end subroutine print_type
end program example_jsonfortran
