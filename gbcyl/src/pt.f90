module pt
  use nrtype
  use mpi
  use particle 
  use class_poly_box
  include 'rng.inc'
  use class_parameter_writer
  use class_parameterizer
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
  public :: pt_writeparameters
  public :: pt_resetcounters
  public :: pt_test_particle_exchange

  interface pt_init
    module procedure pt_initinternal, pt_initwtparameterizer
  end interface

  interface pt_exchange
    module procedure pt_exchange !!, pt_exchange2
  end interface

  integer, save :: dptype
  integer, save :: particletype
  integer, save :: idabove
  integer, save :: idbelow
  integer, save :: id
  integer, save :: ntasks
  integer, save :: idirection = 0
  integer, save :: nacceptedup = 0
  integer, save :: ntrialsup = 0
  integer, save :: confid = -1

contains

!! Creates the MPI types used by this module for communicating configurations 
!! from one process to another. 
!! 
subroutine pt_initinternal()
  integer :: ierr
  integer :: oldtypes(0:1)
  integer :: blockcounts(0:1)
  integer(KIND=MPI_ADDRESS_KIND) :: offsets(0:1)
  integer(KIND=MPI_ADDRESS_KIND) :: extent,lb
  type(particledat) :: aparticle
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr) 
  if(mod(ntasks, 2) /= 0) then
    write(*,*) 'pt_init: This implementation of parallel tempering supports' &
      // ' only an even number of tasks.'
    call MPI_FINALIZE()
    !! This does not stop the program in vuori. Why?
    stop 'Program stopped by pt_init.'
  end if
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr) 
  idabove = mod(id + 1, ntasks)
  idbelow = mod(id - 1 + ntasks, ntasks)
  !! :NOTE: Setting id of particle configuration to id of process. 
  !! :TODO: Make exchange of conf id:s working with restarts.
  if (confid < 0) confid = id
  !if (precision(nrdouble) == precision(somedouble) .and. range(nrdouble) == range(somedouble)) then
  dptype = MPI_REAL8
  !! Create particle type 
  !! Set up description of the coordinate variables.
  call MPI_GET_ADDRESS(aparticle%x, offsets(0), ierr)
  !offsets(0) =  
  oldtypes(0) = MPI_REAL8 
  blockcounts(0) = 6
  !! Setup description of the logical rod variable
  !! Need to first figure offset by getting size of dptype 
  !call MPI_TYPE_GET_EXTENT(dptype, lb, extent, ierr) 
  !if (ierr/=0) stop 'Error: MPI_TYPE_GET_EXTENT failed!!!'
  call MPI_GET_ADDRESS(aparticle%rod, offsets(1), ierr)
  offsets(1)=offsets(1)-offsets(0)
  offsets(0)=0
  !offsets(1) = 6 * extent 
  oldtypes(1) = MPI_LOGICAL8 
  blockcounts(1) = 1 
  !! Now define structured type and commit it  
  call MPI_TYPE_CREATE_STRUCT(2, blockcounts, offsets, oldtypes, particletype, ierr) 
  if (ierr/=0) stop 'Error: Could not create struct for particletype!!!'
  call MPI_TYPE_COMMIT(particletype, ierr) 
  if (ierr/=0) stop 'Error: Could not commit particletype'
  call MPI_TYPE_GET_EXTENT(particletype, lb, extent, ierr)
  !call MPI_TYPE_SIZE(particletype, particlesize, ierr)
  !write(*,*)'Created datatype with extent and size ', extent, sizeof(aparticle)
end subroutine

subroutine setup_test(particles, id)
  type(particledat), intent(inout) :: particles(2)
  integer, intent(in) :: id
  integer :: i
  do i=1,2
    particles(i)%rod=.true.
  end do
  if (mod(id,2)==0) then
    particles(2)%rod=.false.
  end if
end subroutine

subroutine pt_test_particle_exchange 
  integer, parameter :: testtag=77
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr
  integer :: id
  type(particledat) :: particles(2)
  integer, parameter :: n_lots=1000
  type(particledat) :: lots_of_particles(n_lots)
  type(particledat) :: particlebuffer(2)
  integer :: id_other
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)
  if (ierr/=0) write(*, *) 'mpi_comm_rank failed'
  
  if (mod(id,2)==0) then
    id_other=id+1
  else 
    id_other=id-1
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !! Test sending and receiving of particles
  call setup_test(particles, id)

  if (mod(id,2)==0) then
      call MPI_SEND(particles(2), 1, particletype, id_other, testtag, MPI_COMM_WORLD, ierr)
  else
    call MPI_RECV(particles(2), 1, particletype, id_other, testtag, MPI_COMM_WORLD, status, ierr) 
  end if
  if (mod(id,2)==1 .and. particles(2)%rod) then
    write(*,*) 'MPI_SEND(aparticle) + MPI_RECV(aparticle) failed'
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)



  !! Test sending and receiving an array of particles with the regular 
  !! MPI_SEND and MPI_RECV
  call setup_test(particles, id)

  if (mod(id,2)==0) then
      call MPI_SEND(particles, 2, particletype,id_other, testtag, MPI_COMM_WORLD, ierr)
  else
    call MPI_RECV(particles, 2, particletype, id_other, testtag, MPI_COMM_WORLD, status, ierr) 
  end if
  if (mod(id,2)==1 .and. particles(2)%rod) then
    write(*,*) 'MPI_SEND(particles) + MPI_RECV(particles) failed'
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)



  !! Test sending and receiving an array of particles with MPI_SENDRECV
  call setup_test(particles, id)
  particlebuffer=particles
  call MPI_SENDRECV(particlebuffer, 2, particletype,id_other, testtag,&
    & particles, 2, particletype,id_other, testtag, MPI_COMM_WORLD, status,&
    & ierr)
  if (mod(id,2)==1 .and. particles(2)%rod) then
    write(*,*) 'MPI_SENDRECV(particles) failed'
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)



  !! Test sending and receiving an array of particles with MPI_SENDRECV_REPLACE
  call setup_test(particles, id)
  call mpi_sendrecv_replace(particles, 2, particletype, &
      & id_other, testtag, id_other, testtag, MPI_COMM_WORLD, status,&
      & ierr)
  if (ierr/=0) then
    write(*, *) id, 'MPI_SENDRECV_REPLACE(particles) ierr=', ierr
  end if

  if (mod(id,2)==1 .and. particles(2)%rod) then
    write(*,*) 'MPI_SENDRECV_REPLACE(particles) failed'
  !else if(mod(id,2)==1) then
  !  write(*, *) 'MPI_SENDRECV_REPLACE(particles) succeeded'
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !! Test sending and receiving an large array of particles with 
  !! MPI_SENDRECV_REPLACE
  call setup_test(lots_of_particles(1:2), id)
  write(*, *) 'sizeof(lots_of_particles)=',sizeof(lots_of_particles)
  call mpi_sendrecv_replace(lots_of_particles, n_lots, particletype, &
      & id_other, testtag, id_other, testtag, MPI_COMM_WORLD, status,&
      & ierr)
  if (ierr/=0) then
    write(*, *) id, 'MPI_SENDRECV_REPLACE(lots_of_particles) ierr=', ierr
  end if

  if (mod(id,2)==1 .and. particles(2)%rod) then
    write(*,*) 'MPI_SENDRECV_REPLACE(lots_of_particles) failed'
  end if
end subroutine

!! Initializes this module using a parameterizer object to get the parameters
!! for this module. 
!!
!! @p reader the object which reads the parameters. 
!! 
subroutine pt_initwtparameterizer(reader)
  type(parameterizer), intent(in) :: reader
  call getparameter(reader, 'confid', confid)
  call getparameter(reader, 'idirection', idirection)
  call getparameter(reader, 'nacceptedup', nacceptedup)
  call getparameter(reader, 'ntrialsup', ntrialsup)
  call pt_init()
end subroutine

!! Frees the mpi types created by this module. To be called when use of this 
!! module by the program hass ended.
!!
subroutine pt_finalize()
  integer :: ierr
  call mpi_type_free(particletype, ierr)
end subroutine

!! Saves parameters and observables of this module using a parameter_writer 
!! object.
!! 
!! @p pwriter the parameter writer object to be used. 
!! 
subroutine pt_writeparameters(pwriter)
  type(parameter_writer), intent(in) :: pwriter
  call writecomment(pwriter, 'Parallel tempering parameters')
  call writeparameter(pwriter, 'confid', confid)
  call writeparameter(pwriter, 'idirection', idirection)
  call writeparameter(pwriter, 'nacceptedup', nacceptedup)
  call writeparameter(pwriter, 'ntrialsup', ntrialsup)
  call writeparameter(pwriter, 'current_ptaccratioup', &
  real(nacceptedup, dp)/real(ntrialsup, dp))
end subroutine

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
subroutine pt_move(beta, energy, particles, nparticles, simbox, genstate, isaccepted)
  real(dp), intent(in) :: beta
  real(dp), intent(inout) :: energy
  type(particledat), dimension(:), intent(inout) :: particles
  integer, intent(inout) :: nparticles
  type(poly_box), intent(inout) :: simbox
  type(rngstate), intent(inout) :: genstate
  logical, intent(out) :: isaccepted
  real(dp) :: betan
  real(dp) :: energyn
  type(particledat), dimension(nparticles) :: particlesn
  integer :: nparticlesn
  type(poly_box) :: simboxn
  integer :: destid
  real(dp) :: rand 
  real(dp) :: randn
  integer :: confidn
  isaccepted = .false.
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
  rand = rng(genstate)
  !! Make temporary variables. 
  betan = beta
  energyn = energy
  particlesn(1:nparticles) = particles(1:nparticles)
  nparticlesn = nparticles
  simboxn = simbox
  randn = rand
  confidn = confid
  !! Exchange temporary variables
  call pt_exchange(destid, betan, energyn, particlesn, nparticlesn, &
    simboxn, randn, confidn)
  if (mod(id, 2) == 1) rand = randn
  if(exp((beta - betan) * (energy - energyn)) > rand ) then
    energy = energyn
    particles(1:nparticles) = particlesn(1:nparticles)
    nparticles = nparticlesn
    simbox = simboxn
    confid  = confidn
    if (destid == idabove) then
      nacceptedup = nacceptedup + 1
    end if
    isaccepted = .true.
  end if
  if (destid == idabove) then
    ntrialsup = ntrialsup + 1
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
  simbox, rand, confid)
  integer, intent(in) :: destid
  real(dp), intent(inout) :: beta
  real(dp), intent(inout) :: energy
  type(particledat), dimension(:), intent(inout) :: particles
  integer, intent(inout) :: nparticles
  type(poly_box), intent(inout) :: simbox
  real(dp), intent(inout) :: rand
  integer, intent(inout) :: confid
  integer, parameter :: firsttag = 50
  integer, parameter :: secondtag = 52
  integer, dimension(MPI_STATUS_SIZE) :: status 
  integer, parameter :: msgsize = 7
  real(dp), dimension(msgsize) :: msg 
  integer :: ierr
  msg(1) = beta 
  msg(2) = energy 
  msg(3) = rand
  msg(4) = getx(simbox)
  msg(5) = gety(simbox)
  msg(6) = getz(simbox)
  msg(7) = real(confid, dp)
  call mpi_sendrecv_replace(msg(1:msgsize), msgsize, dptype, destid, &
    & firsttag, destid, firsttag, MPI_COMM_WORLD, status, ierr)   
  if (ierr /= 0) then
    stop 'pt: first mpi_sendrecv_replace failed' 
  !else
  !  write(*, *) 'pt: first mpi_sendrecv_replace succeeded'
  end if
  beta = msg(1)
  energy = msg(2)
  rand = msg(3) 
  call setx(simbox, msg(4))
  call sety(simbox, msg(5))
  call setz(simbox, msg(6))
  confid = int(msg(7) + 0.5_dp)

  call mpi_sendrecv_replace(particles, nparticles, &
    & particletype, destid, secondtag, destid, secondtag, MPI_COMM_WORLD, &
    & status, ierr)    
  if (ierr /= 0) then 
    write(*,*) 'pt: second mpi_sendrecv_replace failed with', destid
  !else
  !  write(*, *) 'pt: second mpi_sendrecv_replace succeeded'
  end if
end subroutine

!! Returns the number of accepted pt exchange moves between thisid and
!! thisid + 1.
!! 
pure function pt_getacceptedup() result(naccepted)
  integer :: naccepted 
  naccepted = nacceptedup
end function

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

!! Adjusts the temperatures of the replicas to get evenly large exchange rates.
!! Breaks detailed balance. Do not use during production runs. DEPRECATED.
!! 
!! @p temperatures the array of temperatures to adjust.
!! @p number of accepted swaps between this replica and the replica above it.
!! 
subroutine pt_adjusttemperatures(temperatures, nswaps)
  real(dp), dimension(:), intent(inout) :: temperatures
  integer, dimension(:), intent(in) :: nswaps
  real(dp) :: avgnswaps
  integer :: i
  real(dp), dimension(size(nswaps)) :: shifts
  real(dp) :: shift
  !! Calculate average number of accepted moves
  avgnswaps = real(sum(nswaps(1:ntasks - 1)), dp) / real(ntasks - 1, dp)
  ! Keep endpoint temperatures id == 0 and id == ntasks - 1 fixed.
  do i = 1, ntasks - 1
    shifts(i) = 0.05_dp * (temperatures(i + 1) - temperatures(i))
    if (i > 1 .and. i < ntasks) then
       shift = min(shifts(i), shifts(i - 1))
      if (real(nswaps(i), dp) < avgnswaps) then 
        !! Shift temperature i closer to i + 1
        temperatures(i) = temperatures(i) + shift
      end if
      if (real(nswaps(i - 1), dp) < avgnswaps) then
        !! Shift temperature i closer to i - 1
        temperatures(i) = temperatures(i) - shift
      end if
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
  if (id == 0) then
    allocate(temperatures(ntasks), nswaps(ntasks), stat = astat)
    if (astat /= 0) then 
      stop 'Memory allocation failed in module pt.'
    end if
  end if
  call mpi_gather(pt_getacceptedup(), 1, mpi_integer, nswaps, &
  1, mpi_integer, 0, mpi_comm_world, rc)
  call mpi_gather(temperature, 1, dptype, temperatures, &
  1, dptype, 0, mpi_comm_world, rc)
  !! Make adjustment calculations in id == 0
  if (id == 0) then
    call pt_adjusttemperatures(temperatures, nswaps)
  end if
  call mpi_scatter(temperatures, 1, dptype, temperature,&
  1, dptype, 0, mpi_comm_world, rc) 
end subroutine

subroutine pt_resetcounters
  nacceptedup = 0
  ntrialsup = 0
end subroutine

end module pt
