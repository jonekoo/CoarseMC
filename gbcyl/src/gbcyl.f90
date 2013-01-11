!> This is the main program of the ptgbcyl liquid crystal simulation software.
!! Its purpose is just to call the necessary routines from the simulation 
!! engine which should be exchangeable.
!!
program gbcyl
  use mc_engine, only: init, run, finalize
  use mpi
  use pt
  implicit none
  integer :: ierr
  integer :: id
  integer :: i
  integer :: ntasks
  call mpi_init(ierr)
  if (ierr /= 0) then
    stop 'MPI initialization failed. Check your MPI environment!'
  end if
  call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
  do i=0, ntasks-1
    if (i == id) call init(id, ispt = .true.)
    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end do
  !call pt_test_particle_exchange
  call run 
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call finalize(id)
  call mpi_finalize(ierr)
end program gbcyl
