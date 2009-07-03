module pt
use nrtype
use mpi
use particle 
use box
use mtmod
implicit none



  public :: pt_init
  public :: pt_exchange
  public :: pt_move
  public :: pt_finalize

  private

  integer, save :: dp_type
  integer, save :: particletype
  integer, save :: id_above
  integer, save :: id_below
  integer, save :: id
  integer, save :: ntasks
  integer, save :: i_direction = 0



contains



  subroutine pt_init()
  implicit none
  integer :: ierr
  integer :: oldtypes(0:1)
  integer :: blockcounts(0:1)
  integer :: offsets(0:1)
  integer :: extent
    call MPI_TYPE_CREATE_F90_REAL(dp, MPI_UNDEFINED, dp_type, ierr) 
    !! Set up description of the coordinate variables.
    offsets(0) = 0 
    oldtypes(0) = dp_type 
    blockcounts(0) = 6
   !! Setup description of the logical rod variable
   !! Need to first figure offset by getting size of dp_type 
   call MPI_TYPE_EXTENT(dp_type, extent, ierr) 
   offsets(1) = 6 * extent 
   oldtypes(1) = MPI_LOGICAL 
   blockcounts(1) = 1 
   !! Now define structured type and commit it  
   call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, particletype, ierr) 
   call MPI_TYPE_COMMIT(particletype, ierr) 
   call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr) 
   call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr) 
   id_above = mod(id + 1, ntasks)
   id_below = mod(id - 1 + ntasks, ntasks)
  end subroutine pt_init



  subroutine pt_finalize()
  implicit none
    integer :: ierr
    call mpi_type_free(particletype, ierr)
  end subroutine pt_finalize



  !! Tries to make a parallel tempering Hamiltonian exchange with the process 
  !! in the next lower temperature or the next higher temperature. Direction
  !! is changed every time the routine is being called.   
  !! 
  !! Currently works only if the communicating processes have the same number 
  !! of particles
  !! 
  !! @p beta inverse temperature 1/kT of the calling process. 
  !! @p energy is the total energy (NVT) or enthalpy (NPT) of the process' 
  !!    system.
  !! @p particles array of particles the process is holding. 
  !! @p n_particles dimension of @p particles
  !! @p simbox desribes the simulation cell dimensions and periodicity. 
  !!
  !!
  subroutine pt_move(beta, energy, particles, n_particles, simbox)
  implicit none
  real(dp), intent(in) :: beta
  real(dp), intent(inout) :: energy
  type(particledat), dimension(:), intent(inout) :: particles
  integer, intent(inout) :: n_particles
  type(boxdat), intent(inout) :: simbox
  real(dp) :: beta_n
  real(dp) :: energy_n
  type(particledat), dimension(n_particles) :: particles_n
  integer :: n_particles_n
  type(boxdat) :: simbox_n
  INTEGER :: dest_id
  real(dp) :: acceptance = -1._dp
  real(dp) :: rand 
  real(dp) :: rand_n
    !! Select communication direction
    i_direction = mod(i_direction + 1, 2)
    if (mod(id, 2) == 0) then 
      if(1 == i_direction) then 
        dest_id = id_above
      else 
        dest_id = id_below
      end if
    else 
      if(1 == i_direction) then
        dest_id = id_below
      else 
        dest_id = id_above
      end if
    end if 
    !! Calculate random number for accepting the move.
    rand = grnd()
    !! Make temporary variables. 
    beta_n = beta
    energy_n = energy
    particles_n(1:n_particles) = particles(1:n_particles)
    n_particles_n = n_particles
    simbox_n = simbox
    rand_n = rand
    !! Exchange temporary variables
    call pt_exchange(dest_id, beta_n, energy_n, particles_n, n_particles_n, &
      simbox_n, rand_n)
    if (mod(id, 2) == 1) rand = rand_n
    !! :TODO: Check validity of expression below
    acceptance = exp((beta - beta_n) * (energy - energy_n)) - rand
    if(acceptance > 0._dp) then
      !! Swap:
      energy = energy_n
      particles(1:n_particles) = particles_n(1:n_particles)
      n_particles = n_particles_n
      simbox = simbox_n
    end if
  end subroutine pt_move



  !! Exchanges @p beta, @p energy, @p particles, @p n_particles @p simbox 
  !! dimensions and @p rand with task @p dest_id 
  !!
  !! Currently works only if the communicating processes have the same number 
  !! of particles
  !!
  !! @p dest_id the task id with which the exchange is made with. 
  !! @p beta the inverse temperature. 
  !! @p energy the total energy (NVT) or total enthalpy (NPT) of the process.
  !! @p particles the array of particles. 
  !! @p n_particles the number of particles. 
  !! @p simbox the simulation cell. 
  !! @p rand a random number for accepting the exchange. 
  !! 
  !!
  subroutine pt_exchange(dest_id, beta, energy, particles, n_particles, &
    simbox, rand)
  implicit none
  integer, intent(in) :: dest_id
  real(dp), intent(inout) :: beta
  real(dp), intent(inout) :: energy
  type(particledat), dimension(:), intent(inout) :: particles
  integer, intent(inout) :: n_particles
  type(boxdat), intent(inout) :: simbox
  real(dp), intent(inout) :: rand
  integer, parameter :: first_tag = 50
  integer, parameter :: second_tag = 52
  integer, dimension(MPI_STATUS_SIZE) :: status 
  integer :: msg_size = 6
  real(dp), dimension(6) :: msg 
  integer :: ierr
    msg(1) = beta 
    msg(2) = energy 
    msg(3) = simbox%lx
    msg(4) = simbox%ly
    msg(5) = simbox%lz
    msg(6) = rand
    call mpi_sendrecv_replace(msg(1:msg_size), msg_size, dp_type, dest_id, &
      first_tag, dest_id, first_tag, MPI_COMM_WORLD, status, ierr)
    call mpi_sendrecv_replace(particles(1:n_particles), n_particles, &
      particletype, dest_id, second_tag, dest_id, second_tag, MPI_COMM_WORLD, &
      status, ierr)    
    beta = msg(1)
    energy = msg(2)
    simbox%lx = msg(3)
    simbox%ly = msg(4)
    simbox%lz = msg(5)
    rand = msg(6) 
  end subroutine pt_exchange



  real(dp) function task_temperature(pt_low, pt_high)
  use nrtype
  implicit none
  real(dp), intent(in) :: pt_low
  real(dp), intent(in) :: pt_high
    integer :: n_tasks
    integer :: rc
    integer :: id
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, rc)
    call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc)    
    task_temperature = pt_low + (pt_high - pt_low) * real(id, dp)/ &
      real(n_tasks-1, dp)
  end function task_temperature



end module pt
