
program gbcyl
  use mc_engine, only: init, run, finalize
  use mpi
  implicit none
  integer :: ierr
  call mpi_init(ierr)
  if (ierr /= 0) then
    stop 'MPI initialization failed. Check your MPI environment!'
  end if
  call init
  call run 
  call finalize
  call mpi_finalize(ierr)
end program gbcyl
