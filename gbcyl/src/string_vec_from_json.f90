program string_vec_from_json
  use iso_fortran_env, only: error_unit
  use json_module
  use m_json_wrapper
  implicit none

  character(kind=CK, len=:), allocatable :: vec(:)
  character(kind=CK, len=20), allocatable :: tempvec(:)
  type(json_file) :: jsonf
  type(json_value), pointer :: json_val, jsonvec, child, temp
  logical :: found
  character(kind=CK, len=:), allocatable :: name
  character(len=20) :: justsomeshit
  integer :: i
  !allocate(tempvec(0))
  call json_initialize()
  call jsonf%load_file('test.json')
  call jsonf%get('$', json_val)
  if (.not. associated(json_val)) write(error_unit, *) 'json_val not associated.'
  !call get_string_vec_parameter(json_val, 'participants', vec)
  call json_get(json_val, 'name', name, found)
  if (found) then
     write(*, *) name
  else
     write(error_unit, *) 'name not found'
  end if
  !call json_get(json_val, 'participants', vec, found)
  call json_get(json_val, 'participants', jsonvec)
  if (json_count(jsonvec) < 2) then
     write(error_unit, *) 'less than 2 participants'
  else
     allocate(character(len=20) :: vec(2))
     do i = 1, json_count(jsonvec)
        call json_get_child(jsonvec, i, child)
        call json_get(child, name)
        vec(i) = name
     end do
  end if
  write(*, *) vec
  deallocate(vec)
  allocate(character(kind=CK, len=20) :: vec(0))
  call json_get(json_val, 'participants', temp)
  call json_get(temp, tempvec)
  write(*, *) tempvec
end program string_vec_from_json
