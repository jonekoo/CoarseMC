!> This is the main program of the ptgbcyl simulation software for
!! coarse-grained particles. 
!!
!! @todo create two different simulation programs: with MPI and without
!!       using e.g. conditional compilation.
!!
program gbcyl
  use mc_engine
  use mpi
  use json_module
  use cla
  implicit none
  integer :: ierr
  integer :: id
  integer :: ntasks
  character(len=80) :: input_filename
  character(len=:), allocatable :: parameter_infile
  character(len=80) :: output_filename
  character(len=:), allocatable :: parameter_outfile
  character(len=80) :: restart_filename
  character(len=:), allocatable :: parameter_restartfile
  integer :: nargs
  
  call mpi_init(ierr)
  call json_initialize()
  if (ierr /= 0) then
    stop 'MPI initialization failed. Check your MPI environment!'
  end if
  call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
  !call mce_init(id, ntasks)

  !nargs = command_argument_count()
  !if (nargs < 1) then
  !   write(error_unit, *) 'ERROR: parameter file not given as input.'
  !   write(error_unit, *) 'Examples:'
  !   write(error_unit, *) &
  !        'If the files are input-0000.json input-0001.json etc:'
  !   write(error_unit, *) 'mpirun -np 1 ./ptgbcyl input-_I4_.json'
  !   write(error_unit, *) 'If the files are input-0.json input-1.json etc:'
  !   write(error_unit, *) 'mpirun -np 1 ./ptgbcyl input-_I_.json'
  !   stop
  !end if
  !call get_command_argument(1, input_filename)
  call cla_init
  call cla_register(key='-i', longkey='--input', &
       description='input parameter file', &
       kkind=cla_char, default='input-_I_.json')
  call cla_register(key='-o', longkey='--output', &
       description='output parameter file', &
       kkind=cla_char, default='output-_I_.json')
  call cla_register(key='-r', longkey='--restart-file', &
       description='File to write parameters for restart.', &
       kkind=cla_char, default='restart-_I_.json')
  call cla_validate('ptgbcyl')
  
  call cla_get('--input', input_filename)
  call parse_filename(input_filename, parameter_infile)

  call cla_get('--output', output_filename)
  call parse_filename(output_filename, parameter_outfile)
 
  call cla_get('--restart-file', restart_filename)
  call parse_filename(restart_filename, parameter_restartfile)
 
  call mce_init_json(id, ntasks, parameter_infile, parameter_outfile, &
       parameter_restartfile)

  call run 
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call finalize(id)
  call mpi_finalize(ierr)
contains
  subroutine parse_filename(input_filename, replica_filename)
    character(len=*), intent(in) :: input_filename
    character(len=:), allocatable :: replica_filename
    character(len=20) :: idchar
    integer :: begin
    begin = index(input_filename, '_I_')
    if (begin > 0) then
       write(idchar, '(I20)') id
       replica_filename = input_filename(1:begin - 1) // trim(adjustl(idchar))&
            // input_filename(begin + 3:)
    else
       write(error_unit, *) 'ERROR: Could not convert filename '// &
            input_filename // ' to something meaningful.'
       write(error_unit, *) &
            'Filename should be e.g. input-_I_.json. Here, _I_ is converted ' &
            // ' to be the replica index. So the input is read from ' &
            // 'input-0.json, input-1.json,...'
       stop
    end if
  end subroutine parse_filename
end program gbcyl
