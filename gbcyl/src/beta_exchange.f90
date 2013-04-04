module beta_exchange
use mpi
use nrtype
include 'rng.inc'
implicit none

contains

!> Implements parallel tempering by swapping temperatures.
!! 
!! @param beta the inverse temperature of this replica.
!! @param energy the energy (NVT) or enthalpy (NPT) of this replica.
!! @param p tells n_tasks^p trial swaps are made. n_tasks is the number of 
!! replicas. Chodera et al. recommend p=3-5.
!! @param genstate random number generator state. 
!!
!! @see D. Earl. Parallel tempering: Theory, applications, and new perspectives.
!! Physical Chemistry Chemical Physics, 7(23):3910, 2005. 
!! @see J. D. Chodera and M. R. Shirts. Replica exchange and expanded ensemble 
!! simulations as gibbs sampling: Simple improvements for enhanced mixing. 
!! The Journal of Chemical Physics, 135(19):194110, 2011.
!!
subroutine try_beta_exchanges(beta, energy, p, genstate)
  real(dp), intent(inout) :: beta, energy
  integer, intent(in) :: p
  type(rngstate), intent(inout) :: genstate
  integer :: n_tasks, ierr, task_id, i, j, k
  integer, parameter :: root_id = 0
  real(dp), allocatable :: betas(:), energies(:)
  real(dp) :: r, temp
  call mpi_comm_size(MPI_COMM_WORLD, n_tasks, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, task_id, ierr)
  if (task_id == root_id) then
    allocate(betas(n_tasks), energies(n_tasks))
  end if
  call mpi_gather(beta, 1, MPI_REAL8, betas, 1, MPI_REAL8, root_id, MPI_COMM_WORLD, ierr)
  call mpi_gather(energy, 1, MPI_REAL8, energies, 1, MPI_REAL8, root_id, MPI_COMM_WORLD, ierr)
  if (task_id == root_id) then
    do i = 1, n_tasks**p
      !! Generate random task_ids: 
      call rng(genstate, r)
      j = int(r * real(n_tasks, dp)) + 1
      call rng(genstate, r)
      k = int(r * real(n_tasks, dp)) + 1
      if (j /= k) then
        !! Try to exchange betas:
        call rng(genstate, r)
        if(exp((betas(j) - betas(k)) * (energies(j) - energies(k))) > r ) then
          temp = betas(j)
          betas(j) = betas(k)
          betas(k) = temp
          !write(*, *) 'Made a temperature swap between:', j, k
        end if
      end if
    end do
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_SCATTER(betas, 1, MPI_REAL8, beta, 1, MPI_REAL8, root_id, MPI_COMM_WORLD, ierr)
end subroutine


end module
