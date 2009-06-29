test_suite pt

test pt_wt_rng
use pt
use particle
use box
use nrtype
use mtmod
use mpi
real(dp) :: beta
real(dp) :: energy
type(particledat), dimension(2) :: particles
integer :: n_particles
type(boxdat) :: a_box
integer :: seed = 123456
integer :: ntasks
integer :: id
integer :: rc
logical :: overlap
  !! 1. seed rng
  call sgrnd(seed)
  !! 2. initialize mpi
  call MPI_INIT(rc)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
  if (ntasks /= 2) then
    stop 'This test is suited only for 2 tasks'
  end if
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc) 
  !! 3. initialize configurations
  !! 3.1. Make box
  !! 3.2. task 0: put T-configuration inside the box
  if (id == 0) then
    a_box = make_box(10.0_dp)
    call make_tconf(particles, n_particles)
  else
    a_box = make_box(8.0_dp)
    call make_sidebyside(particles, n_particles) 
  end if
  !! 4. calculate energy
  call pairV(particles(1), particles(2), energy, overlap)
  !! 5. make pt move
  call pt_move(beta, energy, particles, n_particles, a_box)
  !! 6. test if configurations have changed
  !! 7. finalize mpi
  call MPI_FINALIZE(rc)
end test pt_wt_rng

end test_suite pt