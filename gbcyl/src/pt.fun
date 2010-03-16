test_suite pt



test exchange
use pt
use particle
use box
use class_poly_box
use nrtype
use mtmod
use mpi
use gayberne
use configurations
real(dp) :: beta_0 = 0.0_dp
real(dp) :: beta_1 = 1.0_dp
real(dp) :: beta
real(dp) :: energy
type(particledat), dimension(2) :: tconf
type(particledat), dimension(2) :: sidebyside
type(particledat), dimension(2) :: particles
integer :: n_particles
type(poly_box) :: a_box
type(poly_box) :: large_box
type(poly_box) :: small_box
real(dp) :: small_box_side = 10.0_dp
real(dp) :: large_box_side = 20.0_dp
integer :: seed = 123456
integer :: ntasks
integer :: id
integer :: rc
integer :: dest_id
real(dp) :: rand
real(dp) :: rand_0
real(dp) :: rand_1
real(dp) :: E_ss
real(dp) :: E_particles
real(dp) :: E_t
logical :: overlap
  !! Initialize Gay-Berne potential to get the particle configurations right.
  call init(4.4_dp, 20._dp, 1._dp, 1._dp, 1._dp, 1._dp)
  !! 1. seed rng
  call sgrnd(seed)
  rand_0 = grnd()
  rand_1 = grnd()
  !! 2. initialize mpi
  call MPI_INIT(rc)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
  if (ntasks /= 2) then
    stop 'This test is suited only for 2 tasks'
  end if
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc) 
  call pt_init()
  !! 3. initialize configurations
  call make_tconf(tconf, n_particles) 
  !! Make side by side configuration
  call make_sidebyside(sidebyside, n_particles) 
  !! 3.1. Make small box
  small_box = new_box(small_box_side)
  !! Make large box
  large_box = new_box(large_box_side)
  !! 3.2. task 0: put T-configuration inside a small box
  if (id == 0) then
    a_box = small_box
    particles = tconf
    rand = rand_0
    beta = beta_0
  else
  !! Task 1: put side by side configuration to large_box
    a_box = large_box
    particles = sidebyside
    rand = rand_1
    beta = beta_1
  end if
  dest_id = mod(id + 1, 2)
  
  !! Calculate and set potential energy
  call potential((/particles(1)%ux, particles(1)%uy, &
    particles(1)%uz /), (/ particles(2)%ux, particles(2)%uy, &
    particles(2)%uz /), (/particles(1)%x - particles(2)%x, particles(1)%y &
    - particles(2)%y, particles(1)%z - particles(2)%z/), energy, overlap)

  !! 5. make pt exchange
  call pt_exchange(dest_id, beta, energy, particles, n_particles, a_box, rand)
  !! 6. test that configurations have been swapped correctly
  !! Test that particle coordinates have been swapped correctly
  !! Check that box dimensions have been swapped correctly
  !! Check that energies have been swapped correctly
  if (id == 0) then
    !! check that box is now large
    assert_real_equal(large_box_side, get_x(a_box))
    assert_real_equal(large_box_side, get_y(a_box))
    assert_real_equal(large_box_side, get_z(a_box))
    !! We don't test periodicity, cause it is not sent over to the other 
    !! process.
    !! check that configuration is a side by side configuration
    assert_real_equal(particles(1)%x, sidebyside(1)%x)
    assert_real_equal(particles(1)%y, sidebyside(1)%y)
    assert_real_equal(particles(1)%z, sidebyside(1)%z)
    assert_real_equal(particles(1)%ux, sidebyside(1)%ux)
    assert_real_equal(particles(1)%uy, sidebyside(1)%uy)
    assert_real_equal(particles(1)%uz, sidebyside(1)%uz)
    assert_real_equal(particles(2)%x, sidebyside(2)%x)
    assert_real_equal(particles(2)%y, sidebyside(2)%y)
    assert_real_equal(particles(2)%z, sidebyside(2)%z)
    assert_real_equal(particles(2)%ux, sidebyside(2)%ux)
    assert_real_equal(particles(2)%uy, sidebyside(2)%uy)
    assert_real_equal(particles(2)%uz, sidebyside(2)%uz)
    assert_true(particles(1)%rod)
    assert_true(particles(2)%rod)

    !! check that energy is preserved in the transfer
    call potential((/sidebyside(1)%ux, sidebyside(1)%uy, sidebyside(1)%uz &
      /), (/ sidebyside(2)%ux, sidebyside(2)%uy, sidebyside(2)%uz /), &
      (/sidebyside(1)%x - sidebyside(2)%x, sidebyside(1)%y - sidebyside(2)%y,&
       sidebyside(1)%z - sidebyside(2)%z/), E_ss, overlap)
    call potential((/particles(1)%ux, particles(1)%uy, &
      particles(1)%uz /), (/ particles(2)%ux, particles(2)%uy, &
      particles(2)%uz /), (/particles(1)%x - particles(2)%x, particles(1)%y &
      - particles(2)%y, particles(1)%z - particles(2)%z/), E_particles, overlap)
    assert_real_equal(E_ss, E_particles)
    
    !! Check that beta has been exchanged although it is not needed to be 
    !! exchanged.
    assert_real_equal(beta_1, beta)

    !! Check that energy has been exchanged
    assert_real_equal(E_particles, energy)

    !! check that random numbers have been exchanged
    assert_real_equal(rand_1, rand)
    call mpi_barrier(MPI_COMM_WORLD, rc)    !! Makes printing sequential
  else 
    call mpi_barrier(MPI_COMM_WORLD, rc) !! Makes printing sequential

    !! check that box is now small 
    assert_real_equal(small_box_side, get_x(a_box))
    assert_real_equal(small_box_side, get_y(a_box))
    assert_real_equal(small_box_side, get_z(a_box))
    !! We don't test periodicity, cause it is not sent over to the other 
    !! process

    !! check that configuration is now a t-configuration
    assert_real_equal(particles(1)%x, tconf(1)%x)
    assert_real_equal(particles(1)%y, tconf(1)%y)
    assert_real_equal(particles(1)%z, tconf(1)%z)
    assert_real_equal(particles(1)%ux, tconf(1)%ux)
    assert_real_equal(particles(1)%uy, tconf(1)%uy)
    assert_real_equal(particles(1)%uz, tconf(1)%uz)
    assert_true(particles(1)%rod)
    assert_real_equal(particles(2)%x, tconf(2)%x)
    assert_real_equal(particles(2)%y, tconf(2)%y)
    assert_real_equal(particles(2)%z, tconf(2)%z)
    assert_real_equal(particles(2)%ux, tconf(2)%ux)
    assert_real_equal(particles(2)%uy, tconf(2)%uy)
    assert_real_equal(particles(2)%uz, tconf(2)%uz)
    assert_true(particles(2)%rod)

    !! check that energy is preserved in the exchange
    call potential((/tconf(1)%ux, tconf(1)%uy, tconf(1)%uz &
      /), (/ tconf(2)%ux, tconf(2)%uy, tconf(2)%uz /), &
      (/tconf(1)%x - tconf(2)%x, tconf(1)%y - tconf(2)%y,&
       tconf(1)%z - tconf(2)%z/), E_t, overlap)
    call potential((/particles(1)%ux, particles(1)%uy, &
      particles(1)%uz /), (/ particles(2)%ux, particles(2)%uy, &
      particles(2)%uz /), (/particles(1)%x - particles(2)%x, particles(1)%y &
      - particles(2)%y, particles(1)%z - particles(2)%z/), E_particles, overlap)
    assert_real_equal(E_t, E_particles)

    !! Check that energy has been exchanged
    assert_real_equal(E_particles, energy)

    !! Check that beta has been exchanged although it is not needed to be 
    !! exchanged.
    assert_real_equal(beta_0, beta)

    !! Check that random numbers have been exchanged. 
    assert_real_equal(rand_0, rand)
  end if

  call pt_finalize()
  call MPI_FINALIZE(rc)
end test exchange

end test_suite pt