  if (present(error_lb)) then
     if (error_lb > val) then
        write(error_unit, *) 'ERROR: ' // name // ' < ', error_lb, ': ' // &
             name, ' = ', val
        stop
     end if
  end if
  if (present(error_ub)) then
     if (error_ub < val) then
        write(error_unit, *) 'ERROR: ' // name // ' > ', error_ub, ': ' // &
             name, ' = ', val
        stop
     end if
  end if
  if (present(warn_lb)) then
     if (warn_lb > val) then
        write(error_unit, *) 'Warning: ' // name // ' < ', warn_lb, ': ' // &
             name, ' = ', val
     end if
  end if
  if (present(warn_ub)) then
     if (warn_ub < val) then
        write(error_unit, *) 'Warning: ' // name // ' > ', warn_ub, ': ' // &
             name, ' = ', val
     end if
  end if
