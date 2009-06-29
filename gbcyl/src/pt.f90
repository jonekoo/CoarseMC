module pt
use mpi

contains

  subroutine pt_move(beta, energy, particles, n_particles, simbox)
  use nrtype
  use particle
  use box
  implicit none
    real(dp), intent(in) :: beta
    real(dp), intent(inout) :: energy
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(inout) :: n_particles
    type(boxdat), intent(inout) :: simbox
  end subroutine pt_move



  function task_temperature(pt_low, pt_high)
    use nrtype
    implicit none
    real(dp) :: task_temperature
    real(dp), intent(in) :: pt_low
    real(dp), intent(in) :: pt_high
    integer :: n_tasks
    integer :: rc
    integer :: id
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, rc)
    call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc)    
    task_temperature = pt_low + (pt_high - pt_low) * real(id, dp)/real(n_tasks-1, dp)
  end function task_temperature

end module pt
