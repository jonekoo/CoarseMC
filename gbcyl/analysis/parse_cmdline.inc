subroutine parse_cmdline(cmd_line)
  character(len=*), intent(in) :: cmd_line
  character(len=80), allocatable :: args(:)
  character(len=len(cmd_line)) :: message
  integer :: ios
  character(len=len(cmd_line)) :: temp
  !! Remove program name from cmd_line
  call splitstr(cmd_line, " ", args)
  call join(args(2:), "", temp)
  if (trim(adjustl(temp)) == '--help' .or. args(2) == '-h' .or. args(2) == 'help' .or. args(2) == 'h') then
     call print_help
     stop
  else
     temp = "&cmd " // trim(adjustl(temp)) // " /" 
     read(temp, NML = cmd, iostat = ios, iomsg = message)
     if (ios /= 0) then
        write(*, *) "ERROR: ", ios
        write(*, *) message
     end if
  end if
end subroutine



subroutine print_usage
  character(len=80) :: program_name
  call get_command_argument(0, program_name)
  write(*, *) 'Usage: '
  write(*, *) trim(adjustl(program_name)) // ' ' // options_str 
  write(*, *) ''
  write(*, *) 'Call with '
  write(*, *) trim(adjustl(program_name)) // ' --help'
  write(*, *) 'for more information about the options.'
end subroutine print_usage

