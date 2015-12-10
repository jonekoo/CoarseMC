  logical :: status_ok
  character(kind=CK, len=:), allocatable :: error_msg
  call json_check_for_errors(status_ok, error_msg)
  if (.not. status_ok) then
     if (index(error_msg, 'Unable to resolve path') > 0 .or. &
          index(error_msg, 'variable not found') > 0) then
        write(error_unit, *) 'Warning: ' // name // ' not found in ' // &
             json_val%name // '. ' // 'Using default value ' // name //' = ', &
             val
        call json_clear_exceptions()
     else
        write(error_unit, *) 'ERROR: '// error_msg
        stop
     end if
  else
     val = temp
  end if  
