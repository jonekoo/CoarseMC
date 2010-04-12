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
real(dp) :: beta0 = 0.0_dp
real(dp) :: beta1 = 1.0_dp
real(dp) :: beta
real(dp) :: energy
type(particledat), dimension(2) :: tconf
type(particledat), dimension(2) :: sidebyside
type(particledat), dimension(2) :: particles
integer :: nparticles
type(poly_box) :: abox
type(poly_box) :: largebox
type(poly_box) :: smallbox
real(dp) :: smallboxside = 10.0_dp
real(dp) :: largeboxside = 20.0_dp
integer :: seed = 123456
integer :: ntasks
integer :: id
integer :: rc
integer :: destid
real(dp) :: rand
real(dp) :: rand0
real(dp) :: rand1
real(dp) :: Ess
real(dp) :: Eparticles
real(dp) :: Et
logical :: overlap
  !! Initialize Gay-Berne potential to get the particle configurations right.
  call init(4.4_dp, 20._dp, 1._dp, 1._dp, 1._dp, 1._dp)
  !! 1. seed rng
  call sgrnd(seed)
  rand0 = grnd()
  rand1 = grnd()
  !! 2. initialize mpi
  call MPI_INIT(rc)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
  if (ntasks /= 2) then
    stop 'This test is suited only for 2 tasks'
  end if
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc) 
  call pt_init()
  !! 3. initialize configurations
  call maketconf(tconf, nparticles) 
  !! Make side by side configuration
  call makesidebyside(sidebyside, nparticles) 
  !! 3.1. Make small box
  smallbox = new_box(smallboxside)
  !! Make large box
  largebox = new_box(largeboxside)
  !! 3.2. task 0: put T-configuration inside a small box
  if (id == 0) then
    abox = smallbox
    particles = tconf
    rand = rand0
    beta = beta0
  else
  !! Task 1: put side by side configuration to largebox
    abox = largebox
    particles = sidebyside
    rand = rand1
    beta = beta1
  end if
  destid = mod(id + 1, 2)
  
  !! Calculate and set potential energy
  call potential((/particles(1)%ux, particles(1)%uy, &
    particles(1)%uz /), (/ particles(2)%ux, particles(2)%uy, &
    particles(2)%uz /), (/particles(1)%x - particles(2)%x, particles(1)%y &
    - particles(2)%y, particles(1)%z - particles(2)%z/), energy, overlap)

  !! 5. make pt exchange
  call pt_exchange(destid, beta, energy, particles, nparticles, abox, rand)
  !! 6. test that configurations have been swapped correctly
  !! Test that particle coordinates have been swapped correctly
  !! Check that box dimensions have been swapped correctly
  !! Check that energies have been swapped correctly
  if (id == 0) then
    !! check that box is now large
    assert_real_equal(largeboxside, getx(abox))
    assert_real_equal(largeboxside, gety(abox))
    assert_real_equal(largeboxside, getz(abox))
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
       sidebyside(1)%z - sidebyside(2)%z/), Ess, overlap)
    call potential((/particles(1)%ux, particles(1)%uy, &
      particles(1)%uz /), (/ particles(2)%ux, particles(2)%uy, &
      particles(2)%uz /), (/particles(1)%x - particles(2)%x, particles(1)%y &
      - particles(2)%y, particles(1)%z - particles(2)%z/), Eparticles, overlap)
    assert_real_equal(Ess, Eparticles)
    
    !! Check that beta has been exchanged although it is not needed to be 
    !! exchanged.
    assert_real_equal(beta1, beta)

    !! Check that energy has been exchanged
    assert_real_equal(Eparticles, energy)

    !! check that random numbers have been exchanged
    assert_real_equal(rand1, rand)
    call mpi_barrier(MPI_COMM_WORLD, rc)    !! Makes printing sequential
  else 
    call mpi_barrier(MPI_COMM_WORLD, rc) !! Makes printing sequential

    !! check that box is now small 
    assert_real_equal(smallboxside, getx(abox))
    assert_real_equal(smallboxside, gety(abox))
    assert_real_equal(smallboxside, getz(abox))
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
       tconf(1)%z - tconf(2)%z/), Et, overlap)
    call potential((/particles(1)%ux, particles(1)%uy, &
      particles(1)%uz /), (/ particles(2)%ux, particles(2)%uy, &
      particles(2)%uz /), (/particles(1)%x - particles(2)%x, particles(1)%y &
      - particles(2)%y, particles(1)%z - particles(2)%z/), Eparticles, overlap)
    assert_real_equal(Et, Eparticles)

    !! Check that energy has been exchanged
    assert_real_equal(Eparticles, energy)

    !! Check that beta has been exchanged although it is not needed to be 
    !! exchanged.
    assert_real_equal(beta0, beta)

    !! Check that random numbers have been exchanged. 
    assert_real_equal(rand0, rand)
  end if

  call pt_finalize()
  call MPI_FINALIZE(rc)
end test exchange

end test_suite pt