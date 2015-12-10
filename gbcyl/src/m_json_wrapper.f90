module m_json_wrapper
  use iso_fortran_env, only: error_unit
  use json_module
  use num_kind, only: dp
  implicit none
  
  interface get_parameter
     module procedure get_integer_parameter, get_real_parameter, &
          get_logical_parameter, get_string_vec_parameter
  end interface get_parameter

  interface process_not_found
     module procedure process_integer_not_found, process_real_not_found, &
          process_logical_not_found, process_string_vec_not_found
  end interface process_not_found
  
contains

  subroutine get_integer_parameter(json_val, name, val, error_lb, error_ub, &
       warn_lb, warn_ub)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    integer, intent(inout) :: val
    integer, intent(in), optional :: error_lb, error_ub, warn_lb, warn_ub
    integer :: temp
    call json_get(json_val, name, temp)
    call process_not_found(json_val, name, val, temp)
    include 'check_parameter.inc'
  end subroutine get_integer_parameter

  subroutine process_integer_not_found(json_val, name, val, temp)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    integer, intent(inout) :: val
    integer, intent(in) :: temp
    include 'process_not_found.inc'    
  end subroutine process_integer_not_found
  
  subroutine get_real_parameter(json_val, name, val, error_lb, error_ub, &
       warn_lb, warn_ub)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    real(dp), intent(inout) :: val
    real(dp), intent(in), optional :: error_lb, error_ub, warn_lb, warn_ub
    real(dp) :: temp
    call json_get(json_val, name, temp)
    call process_not_found(json_val, name, val, temp)
    include 'check_parameter.inc'
  end subroutine get_real_parameter

  subroutine process_real_not_found(json_val, name, val, temp)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    real(dp), intent(inout) :: val
    real(dp), intent(in) :: temp
    include 'process_not_found.inc'    
  end subroutine process_real_not_found
  
  subroutine get_logical_parameter(json_val, name, val)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    logical, intent(inout) :: val
    logical :: temp
    call json_get(json_val, name, temp)
    call process_not_found(json_val, name, val, temp)
  end subroutine get_logical_parameter
  
  subroutine process_logical_not_found(json_val, name, val, temp)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    logical, intent(inout) :: val
    logical, intent(in) :: temp
    include 'process_not_found.inc'    
  end subroutine process_logical_not_found

  subroutine get_string_vec_parameter(json_val, name, val)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    character(kind=CK, len=*), allocatable, intent(inout) :: val(:)
    type(json_value), pointer :: temp
    character(kind=CK, len=len(val)), allocatable :: vec(:)
    call json_get(json_val, name, temp)
    call json_get(temp, vec)
    call process_not_found(json_val, name, val, vec)
    !call json_destroy(temp, .false.)
  end subroutine get_string_vec_parameter

  subroutine process_string_vec_not_found(json_val, name, val, temp)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    character(kind=CK, len=*), allocatable, intent(inout) :: val(:)
    character(kind=CK, len=*), intent(in) :: temp(:)
    include 'process_not_found.inc'
  end subroutine process_string_vec_not_found
  
end module m_json_wrapper
