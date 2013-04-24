!> This is the main program of the ptgbcyl liquid crystal simulation software.
!! Its purpose is just to call the necessary routines from the simulation 
!! engine which should be exchangeable.
!!
program gbcyl
  use mc_engine, only: init, run, finalize
  use mpi
  implicit none
  integer :: ierr
  integer :: id
  integer :: ntasks
  call mpi_init(ierr)
  if (ierr /= 0) then
    stop 'MPI initialization failed. Check your MPI environment!'
  end if
  call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
  !! call init(input_configuration_reader, input_parameters_reader)
  !do i=0, ntasks-1
  call init(id, ntasks)
  !end do
  !call pt_test_particle_exchange
  ! call sample(n_equilibration, output_configuration_writer, output_parameters_writer, restart_configuration_writer, restart_parameters_writer, adjust=.true.)
  ! call sample(n_production, output_configuration_writer, output_parameters_writer, restart_configuration_writer, restart_parameters_writer)
  call run 
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call finalize(id)
  call mpi_finalize(ierr)
end program gbcyl
