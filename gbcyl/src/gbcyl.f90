!> This is the main program of the ptgbcyl simulation software for
!! anisotropic particles. 
!!
!! @todo create two different simulation programs: with MPI and without
!!       using e.g. conditional compilation.
!!
program gbcyl
  use mc_engine
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
  call mce_init(id, ntasks)
  call run 
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call finalize(id)
  call mpi_finalize(ierr)
end program gbcyl
