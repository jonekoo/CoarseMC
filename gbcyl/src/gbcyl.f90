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
  call init(id, ntasks)
  call run 
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call finalize(id)
  call mpi_finalize(ierr)
end program gbcyl
