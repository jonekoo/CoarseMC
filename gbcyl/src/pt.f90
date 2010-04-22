module pt
  use nrtype
  use mpi
  use particle 
  use class_poly_box
  use mtmod
  implicit none
  private

  public :: pt_init
  public :: pt_exchange
  public :: pt_move
  public :: pt_finalize
  public :: pt_temperature
  public :: pt_getacceptedup
  public :: pt_setacceptedup
  public :: pt_adjusttemperature

  integer, save :: dptype
  integer, save :: particletype
  integer, save :: idabove
  integer, save :: idbelow
  integer, save :: id
  integer, save :: ntasks
  integer, save :: idirection = 0
  integer, save :: nacceptedup = 0

contains

  subroutine pt_init()
    integer :: ierr
    integer :: oldtypes(0:1)
    integer :: blockcounts(0:1)
    integer :: offsets(0:1)
    integer :: extent
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr) 
    if(mod(ntasks, 2) /= 0) then
      write(*,*) 'pt_init: parallel tempering supports only an even number ', &
        'of tasks.'
      call MPI_FINALIZE()
      stop 'Program stopped by pt_init.'
    end if
    call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr) 
    idabove = mod(id + 1, ntasks)
    idbelow = mod(id - 1 + ntasks, ntasks)
    call MPI_TYPE_CREATE_F90_REAL(dp, MPI_UNDEFINED, dptype, ierr) 
  
    !! Create particle type 
    !! Set up description of the coordinate variables.
    offsets(0) = 0 
    oldtypes(0) = dptype 
    blockcounts(0) = 6
    !! Setup description of the logical rod variable
    !! Need to first figure offset by getting size of dptype 
    call MPI_TYPE_EXTENT(dptype, extent, ierr) 
    offsets(1) = 6 * extent 
    oldtypes(1) = MPI_LOGICAL 
    blockcounts(1) = 1 
    !! Now define structured type and commit it  
    call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, particletype, ierr) 
    call MPI_TYPE_COMMIT(particletype, ierr) 
  end subroutine pt_init



  subroutine pt_finalize()
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
  !! @p nparticles dimension of @p particles
  !! @p simbox desribes the simulation cell dimensions and periodicity. 
  !!
  !! :TODO: Make pt_move independent of MPI and only pt_exchange dependent
  !!
  subroutine pt_move(beta, energy, particles, nparticles, simbox)
    real(dp), intent(in) :: beta
    real(dp), intent(inout) :: energy
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(inout) :: nparticles
    type(poly_box), intent(inout) :: simbox
    real(dp) :: betan
    real(dp) :: energyn
    type(particledat), dimension(nparticles) :: particlesn
    integer :: nparticlesn
    type(poly_box) :: simboxn
    integer :: destid
    real(dp) :: rand 
    real(dp) :: randn
    !! Select communication direction
    idirection = mod(idirection + 1, 2)
    if (mod(id, 2) == 0) then 
      if(1 == idirection) then 
        destid = idabove
      else 
        destid = idbelow
      end if
    else 
      if(1 == idirection) then
        destid = idbelow
      else 
        destid = idabove
      end if
    end if 
    !! Calculate random number for accepting the move.
    rand = grnd()
    !! Make temporary variables. 
    betan = beta
    energyn = energy
    particlesn(1:nparticles) = particles(1:nparticles)
    nparticlesn = nparticles
    simboxn = simbox
    randn = rand
    !! Exchange temporary variables
    call pt_exchange(destid, betan, energyn, particlesn, nparticlesn, &
      simboxn, randn)
    if (mod(id, 2) == 1) rand = randn
    !! :TODO: Check validity of expression below
    if(exp((beta - betan) * (energy - energyn)) > rand ) then
      energy = energyn
      particles(1:nparticles) = particlesn(1:nparticles)
      nparticles = nparticlesn
      simbox = simboxn
      if (destid == idabove) then
        nacceptedup = nacceptedup + 1
      end if
    end if
  end subroutine


  !! Exchanges @p beta, @p energy, @p particles, @p nparticles @p simbox 
  !! dimensions and @p rand with task @p destid 
  !!
  !! Currently works only if the communicating processes have the same number 
  !! of particles
  !!
  !! @p destid the task id with which the exchange is made with. 
  !! @p beta the inverse temperature. 
  !! @p energy the total energy (NVT) or total enthalpy (NPT) of the process.
  !! @p particles the array of particles. 
  !! @p nparticles the number of particles. 
  !! @p simbox the simulation cell. 
  !! @p rand a random number for accepting the exchange. 
  !! 
  !! :TODO: Test that the routine is capable of sending and receiving a 
  !! :TODO: different amount of particles. If needed use mpi_get_count 
  !! :TODO: to implement this feature.
  !!
  subroutine pt_exchange(destid, beta, energy, particles, nparticles, &
    simbox, rand)
    integer, intent(in) :: destid
    real(dp), intent(inout) :: beta
    real(dp), intent(inout) :: energy
    type(particledat), dimension(:), intent(inout) :: particles
    integer, intent(inout) :: nparticles
    type(poly_box), intent(inout) :: simbox
    real(dp), intent(inout) :: rand
    integer, parameter :: firsttag = 50
    integer, parameter :: secondtag = 52
    integer, dimension(MPI_STATUS_SIZE) :: status 
    integer :: msgsize = 6
    real(dp), dimension(6) :: msg 
    integer :: ierr
    msg(1) = beta 
    msg(2) = energy 
    msg(3) = rand
    msg(4) = getx(simbox)
    msg(5) = gety(simbox)
    msg(6) = getz(simbox)
    call mpi_sendrecv_replace(msg(1:msgsize), msgsize, dptype, destid, &
      firsttag, destid, firsttag, MPI_COMM_WORLD, status, ierr)
    call mpi_sendrecv_replace(particles(1:nparticles), nparticles, &
      particletype, destid, secondtag, destid, secondtag, MPI_COMM_WORLD, &
      status, ierr)    
    beta = msg(1)
    energy = msg(2)
    rand = msg(3) 
    call setx(simbox, msg(4))
    call sety(simbox, msg(5))
    call setz(simbox, msg(6))
  end subroutine pt_exchange

  !! Returns the number of accepted pt exchange moves between thisid and
  !! thisid + 1.
  !! 
  pure function pt_getacceptedup() result(naccepted)
    integer :: naccepted 
    naccepted = nacceptedup
  end function

  !! Resets counter nacceptedup to zero.
  !!
  subroutine pt_setacceptedup(naccepted)
    integer, intent(in) :: naccepted
    nacceptedup = naccepted
  end subroutine

  real(dp) function pt_temperature(ptlow, pthigh)
    real(dp), intent(in) :: ptlow
    real(dp), intent(in) :: pthigh
    integer :: ntasks
    integer :: rc
    integer :: id
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
    call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc)    
    pt_temperature = ptlow + (pthigh - ptlow) * real(id, dp)/ &
    real(ntasks - 1, dp)
  end function

subroutine pt_adjusttemperatures(temperatures, nswaps)
  real(dp), dimension(:), intent(inout) :: temperatures
  integer, dimension(:), intent(in) :: nswaps
  real(dp) :: avgnswaps
  integer :: i
  real(dp) :: shift
    !! Calculate average number of accepted moves
    avgnswaps = real(sum(nswaps(1:ntasks - 1)), dp) / real(ntasks - 1, dp)   
    !! Keep endpoint temperatures id == 0 and id == ntasks - 1 fixed.
    do i = 1, ntasks - 1
      shift = 0.1_dp * (temperatures(i + 1) - temperatures(i))
      !! Shift temperatures further apart from each other or
      if (real(nswaps(i), dp) < avgnswaps) then 
        !! Shift temperatures closer to each other.
        shift = -shift
      end if
      if (i > 1) then
        temperatures(i) = temperatures(i) - shift
      end if
      if (i + 1 < ntasks) then
        temperatures(i + 1) = temperatures(i + 1) + shift
      end if
    end do
end subroutine

subroutine pt_adjusttemperature(temperature)
!! Adjusting constant could be for example 1/10th of the interval between 
!! initial or instantaneous temperatures. Using interval of instantaneous 
!! temperatures would ensure that the temperatures never cross each other or 
!! overlap. The constant could also be based on distance between averages of 
!! enthalpy, but that would complicate things a lot. 
!! 
  real(dp), intent(inout) :: temperature
  real(dp), dimension(:), allocatable :: temperatures
  integer, dimension(:), allocatable :: nswaps
  integer :: rc
  integer :: astat
!  call mpi_comm_rank(mpi_comm_world, id, rc)
!  call mpi_comm_rank(mpi_comm_world, ntasks, rc)
  if (id == 0) then
    allocate(temperatures(0:ntasks - 1), nswaps(0:ntasks - 1), stat = astat)
    if (astat /= 0) then 
      stop 'Memory allocation failed in module pt.'
    end if
  end if
  call mpi_gather(pt_getacceptedup(), 1, mpi_integer, nswaps(0:ntasks - 1), &
  1, mpi_integer, 0, mpi_comm_world, rc)
  call mpi_gather(temperature, 1, dptype, temperatures(0:ntasks - 1), &
  1, dptype, 0, mpi_comm_world, rc)
  !! Make adjustment calculations in id == 0
  if (id == 0) then
    call pt_adjusttemperatures(temperatures(0:ntasks - 1), &
    nswaps(0:ntasks - 1))
  end if
  call mpi_scatter(temperatures(0:ntasks - 1), 1, dptype, temperature,&
  1, dptype, 0, mpi_comm_world, rc) 
  call pt_setacceptedup(0)
end subroutine

end module pt
