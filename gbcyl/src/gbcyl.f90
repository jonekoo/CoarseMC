program gbcyl
  use mc_engine, only: init, run, finalize, equilibration_sweeps, &
    & production_sweeps, read_restart
  implicit none
  intrinsic command_argument_count, get_command_argument
  character(len=32) :: argument
  character(len=32) :: start_type
  integer :: n_equilibration_sweeps
  integer :: n_production_sweeps
  integer :: n_args
  n_args = command_argument_count()
  if(n_args > 0) then
    call get_command_argument(1, argument) 
    start_type = argument
    if(start_type == "restart") then
      call read_restart   
      if(n_args == 3) then
        call get_command_argument(2, argument)
        read(argument, *) n_equilibration_sweeps
        call equilibration_sweeps(n_equilibration_sweeps)
        call get_command_argument(3, argument)
        read(argument, *) n_production_sweeps
        call production_sweeps(n_production_sweeps)
      end if
    else
      !! Print usage
    end if
  else 
    call init
  end if
  call run 
  call finalize
end program gbcyl


