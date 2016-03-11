!> Implements the parallel tempering algorithm via exchange of inverse
!! temperatures (betas) with MPI.
!!
module beta_exchange
use mpi
use num_kind
use utils
use mt_stream, only: rngstate => mt_state, rng => genrand_double1_s
implicit none

public :: init, finalize, be_reset_counters, try_beta_exchanges, write_stats, &
     be_get_acceptance_ratios

private

integer, allocatable, save :: temperature_index(:)
integer, allocatable, save :: n_accepted(:,:)
integer, allocatable, save :: n_trials(:,:)
integer, parameter :: root_id = 0

contains


!> Initializes the module variables.
!!
!! @pre MPI must be initialized.
!!
!! @param beta the inverse temperature of the MPI task calling this.
!!
subroutine init(beta)
  real(dp), intent(in) :: beta
  integer :: n_tasks, task_id, ierr
  real(dp), allocatable :: betas(:)
  integer :: i, j, data
  real(dp) :: key
  integer, allocatable :: indices(:)

  call mpi_comm_size(MPI_COMM_WORLD, n_tasks, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, task_id, ierr)
  if (task_id == root_id) then
    allocate(n_accepted(n_tasks, n_tasks), &
      n_trials(n_tasks, n_tasks), temperature_index(n_tasks))
    call be_reset_counters()
  end if

  !! Find out the original indices of temperatures
  !! First gather all temperatures to root task
  call mpi_comm_size(MPI_COMM_WORLD, n_tasks, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, task_id, ierr)
  if (task_id == root_id) then
    allocate(betas(n_tasks), indices(n_tasks))
  end if
  call mpi_gather(beta, 1, MPI_REAL8, betas, 1, MPI_REAL8, root_id, &
    MPI_COMM_WORLD, ierr)
  if (task_id == root_id) then
    !! Determine the replica index - temperature index mapping.
    !! First sort replica indices by temperature.
    !!
    !! Insertion sort from T. Cormen et al. here except sort to decreasing order.
    indices(1) = 1
    do j = 2, n_tasks
      key = betas(j)
      data = j
      i = j - 1
      do while(i > 0 .and. betas(i) < key)
        betas(i + 1) = betas(i)
        indices(i + 1) = indices(i)
        i = i - 1
      end do
      betas(i + 1) = key
      indices(i + 1) = data
    end do
    !! After the sorting, indices contains indices of the replicas ordered by
    !! temperature. We need indices of temperatures ordered by replica, so
    do i = 1, n_tasks
      temperature_index(indices(i)) = i
    end do
  end if  
  if (allocated(indices)) deallocate(indices)
  if (allocated(betas)) deallocate(betas)
end subroutine


!> Deallocates the module variables.
subroutine finalize
  if (allocated(temperature_index)) deallocate(temperature_index)
  if (allocated(n_trials)) deallocate(n_trials)
  if (allocated(n_accepted)) deallocate(n_accepted)
end subroutine

!> Resets the counters for temperature swap trials and accepted swaps.
subroutine be_reset_counters()
  n_accepted = 0
  n_trials = 0
end subroutine

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
  real(dp), intent(inout) :: beta
  real(dp), intent(in) :: energy
  integer, intent(in) :: p
  type(rngstate), intent(inout) :: genstate
  integer :: n_tasks, ierr, task_id, i, j, k
  real(dp), allocatable :: betas(:), energies(:)
  real(dp) :: r, temp
  integer :: temp_index
  call mpi_comm_size(MPI_COMM_WORLD, n_tasks, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, task_id, ierr)
  if (task_id == root_id) then
    allocate(betas(n_tasks), energies(n_tasks))
  end if
  call mpi_gather(beta, 1, MPI_REAL8, betas, 1, MPI_REAL8, root_id, &
    MPI_COMM_WORLD, ierr)
  call mpi_gather(energy, 1, MPI_REAL8, energies, 1, MPI_REAL8, root_id, &
    MPI_COMM_WORLD, ierr)
  if (task_id == root_id) then
    do i = 1, n_tasks**p
      !! Generate random task_ids: 
      call rng(genstate, r)
      j = int(r * real(n_tasks, dp)) + 1
      call rng(genstate, r)
      k = int(r * real(n_tasks, dp)) + 1
      if (j /= k) then
        !! Try to exchange betas:
        n_trials(temperature_index(j), temperature_index(k)) = &
          n_trials(temperature_index(j), temperature_index(k)) + 1
        call rng(genstate, r)
        if(exp((betas(j) - betas(k)) * (energies(j) - energies(k))) > r ) then
          temp = betas(j)
          betas(j) = betas(k)
          betas(k) = temp
          !! Book keeping for stats:
          n_accepted(temperature_index(j), temperature_index(k)) = &
            n_accepted(temperature_index(j), temperature_index(k)) + 1
          !! Swap indices in temperature_index
          temp_index = temperature_index(k)
          temperature_index(k) = temperature_index(j)
          temperature_index(j) = temp_index
        end if
      end if
    end do
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_SCATTER(betas, 1, MPI_REAL8, beta, 1, MPI_REAL8, root_id, &
    MPI_COMM_WORLD, ierr)
  if (allocated(betas)) deallocate(betas)
  if (allocated(energies)) deallocate(energies)
end subroutine


!> Writes simple statistics about the acceptance of swaps to the output
!! @p unit.
!!
!! @param unit the output unit to write the statistics to.
!!
subroutine write_stats(unit)
  integer, intent(in) :: unit
  integer :: j, k
  real(dp) :: ratio
  do j = 1, size(n_trials, 1)
    do k = 1, size(n_trials, 2)
      if (n_trials(j, k) > 0) then
        ratio = real(n_accepted(j, k), dp) / real(n_trials(j, k), dp)
      else 
        ratio = -1._dp
      end if
      write(unit, ADVANCE='NO', fmt='(' // fmt_char_dp() // ',1X)') ratio
    end do
    write(unit, *) ''
  end do
end subroutine


!> Returns acceptance ratios of temperature swaps for each replica with
!! all the other replicas.
subroutine be_get_acceptance_ratios(ratios)
  real(dp), allocatable, intent(out) :: ratios(:, :)
  ratios = real(n_accepted, dp)
  where(n_trials > 0)
     ratios = ratios / n_trials
  elsewhere
     ratios = -1._dp
  end where
end subroutine

end module
