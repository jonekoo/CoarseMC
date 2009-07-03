! pt_fun.f90 - a unit test suite for pt.f90
!
! funit generated this file from pt.fun
! at Fri Jul 03 16:13:49 +0300 2009

module pt_fun

 use pt

 implicit none

 logical :: noAssertFailed

 public :: test_pt

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0





 contains

 subroutine exchange

use pt
use particle
use box
use nrtype
use mtmod
use mpi
use gayberne
real(dp) :: beta_0 = 0.0_dp
real(dp) :: beta_1 = 1.0_dp
real(dp) :: beta
real(dp) :: energy
type(particledat), dimension(2) :: tconf
type(particledat), dimension(2) :: sidebyside
type(particledat), dimension(2) :: particles
integer :: n_particles
type(boxdat) :: a_box
type(boxdat) :: large_box
type(boxdat) :: small_box
real(dp) :: small_box_side = 10.0_dp
real(dp) :: large_box_side = 20.0_dp
integer :: seed = 123456
integer :: ntasks
integer :: id
integer :: rc
logical :: overlap
integer :: dest_id
real(dp) :: rand
real(dp) :: rand_0
real(dp) :: rand_1
real(dp) :: E_ss
real(dp) :: E_particles
real(dp) :: E_t
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
  small_box = make_box(small_box_side)
  !! Make large box
  large_box = make_box(large_box_side)
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
  E_particles = potential((/particles(1)%ux, particles(1)%uy, &
    particles(1)%uz /), (/ particles(2)%ux, particles(2)%uy, &
    particles(2)%uz /), (/particles(1)%x - particles(2)%x, particles(1)%y &
    - particles(2)%y, particles(1)%z - particles(2)%z/))
  energy = E_particles 

  !! 5. make pt exchange
  call pt_exchange(dest_id, beta, energy, particles, n_particles, a_box, rand)
  !! 6. test that configurations have been swapped correctly
  !! Test that particle coordinates have been swapped correctly
  !! Check that box dimensions have been swapped correctly
  !! Check that energies have been swapped correctly
  if (id == 0) then
    !! check that box is now large
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (large_box_side &
        +2*spacing(real(large_box_side)) ) &
        .ge. &
        (box_x(a_box)) &
            .and. &
     (large_box_side &
      -2*spacing(real(large_box_side)) ) &
      .le. &
       (box_x(a_box)) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:90]"
      print *, "  ", "box_x(a_box) (", &
 box_x(a_box), &
  ") is not", &
 large_box_side,&
 "within", &
  2*spacing(real(large_box_side))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (large_box_side &
        +2*spacing(real(large_box_side)) ) &
        .ge. &
        (box_y(a_box)) &
            .and. &
     (large_box_side &
      -2*spacing(real(large_box_side)) ) &
      .le. &
       (box_y(a_box)) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:91]"
      print *, "  ", "box_y(a_box) (", &
 box_y(a_box), &
  ") is not", &
 large_box_side,&
 "within", &
  2*spacing(real(large_box_side))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (large_box_side &
        +2*spacing(real(large_box_side)) ) &
        .ge. &
        (box_z(a_box)) &
            .and. &
     (large_box_side &
      -2*spacing(real(large_box_side)) ) &
      .le. &
       (box_z(a_box)) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:92]"
      print *, "  ", "box_z(a_box) (", &
 box_z(a_box), &
  ") is not", &
 large_box_side,&
 "within", &
  2*spacing(real(large_box_side))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    !! We don't test periodicity, cause it is not sent over to the other 
    !! process.
    !! check that configuration is a side by side configuration
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%x &
        +2*spacing(real(particles(1)%x)) ) &
        .ge. &
        (sidebyside(1)%x) &
            .and. &
     (particles(1)%x &
      -2*spacing(real(particles(1)%x)) ) &
      .le. &
       (sidebyside(1)%x) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:96]"
      print *, "  ", "sidebyside(1)%x (", &
 sidebyside(1)%x, &
  ") is not", &
 particles(1)%x,&
 "within", &
  2*spacing(real(particles(1)%x))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%y &
        +2*spacing(real(particles(1)%y)) ) &
        .ge. &
        (sidebyside(1)%y) &
            .and. &
     (particles(1)%y &
      -2*spacing(real(particles(1)%y)) ) &
      .le. &
       (sidebyside(1)%y) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:97]"
      print *, "  ", "sidebyside(1)%y (", &
 sidebyside(1)%y, &
  ") is not", &
 particles(1)%y,&
 "within", &
  2*spacing(real(particles(1)%y))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%z &
        +2*spacing(real(particles(1)%z)) ) &
        .ge. &
        (sidebyside(1)%z) &
            .and. &
     (particles(1)%z &
      -2*spacing(real(particles(1)%z)) ) &
      .le. &
       (sidebyside(1)%z) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:98]"
      print *, "  ", "sidebyside(1)%z (", &
 sidebyside(1)%z, &
  ") is not", &
 particles(1)%z,&
 "within", &
  2*spacing(real(particles(1)%z))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%ux &
        +2*spacing(real(particles(1)%ux)) ) &
        .ge. &
        (sidebyside(1)%ux) &
            .and. &
     (particles(1)%ux &
      -2*spacing(real(particles(1)%ux)) ) &
      .le. &
       (sidebyside(1)%ux) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:99]"
      print *, "  ", "sidebyside(1)%ux (", &
 sidebyside(1)%ux, &
  ") is not", &
 particles(1)%ux,&
 "within", &
  2*spacing(real(particles(1)%ux))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%uy &
        +2*spacing(real(particles(1)%uy)) ) &
        .ge. &
        (sidebyside(1)%uy) &
            .and. &
     (particles(1)%uy &
      -2*spacing(real(particles(1)%uy)) ) &
      .le. &
       (sidebyside(1)%uy) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:100]"
      print *, "  ", "sidebyside(1)%uy (", &
 sidebyside(1)%uy, &
  ") is not", &
 particles(1)%uy,&
 "within", &
  2*spacing(real(particles(1)%uy))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%uz &
        +2*spacing(real(particles(1)%uz)) ) &
        .ge. &
        (sidebyside(1)%uz) &
            .and. &
     (particles(1)%uz &
      -2*spacing(real(particles(1)%uz)) ) &
      .le. &
       (sidebyside(1)%uz) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:101]"
      print *, "  ", "sidebyside(1)%uz (", &
 sidebyside(1)%uz, &
  ") is not", &
 particles(1)%uz,&
 "within", &
  2*spacing(real(particles(1)%uz))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%x &
        +2*spacing(real(particles(2)%x)) ) &
        .ge. &
        (sidebyside(2)%x) &
            .and. &
     (particles(2)%x &
      -2*spacing(real(particles(2)%x)) ) &
      .le. &
       (sidebyside(2)%x) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:102]"
      print *, "  ", "sidebyside(2)%x (", &
 sidebyside(2)%x, &
  ") is not", &
 particles(2)%x,&
 "within", &
  2*spacing(real(particles(2)%x))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%y &
        +2*spacing(real(particles(2)%y)) ) &
        .ge. &
        (sidebyside(2)%y) &
            .and. &
     (particles(2)%y &
      -2*spacing(real(particles(2)%y)) ) &
      .le. &
       (sidebyside(2)%y) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:103]"
      print *, "  ", "sidebyside(2)%y (", &
 sidebyside(2)%y, &
  ") is not", &
 particles(2)%y,&
 "within", &
  2*spacing(real(particles(2)%y))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%z &
        +2*spacing(real(particles(2)%z)) ) &
        .ge. &
        (sidebyside(2)%z) &
            .and. &
     (particles(2)%z &
      -2*spacing(real(particles(2)%z)) ) &
      .le. &
       (sidebyside(2)%z) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:104]"
      print *, "  ", "sidebyside(2)%z (", &
 sidebyside(2)%z, &
  ") is not", &
 particles(2)%z,&
 "within", &
  2*spacing(real(particles(2)%z))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%ux &
        +2*spacing(real(particles(2)%ux)) ) &
        .ge. &
        (sidebyside(2)%ux) &
            .and. &
     (particles(2)%ux &
      -2*spacing(real(particles(2)%ux)) ) &
      .le. &
       (sidebyside(2)%ux) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:105]"
      print *, "  ", "sidebyside(2)%ux (", &
 sidebyside(2)%ux, &
  ") is not", &
 particles(2)%ux,&
 "within", &
  2*spacing(real(particles(2)%ux))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%uy &
        +2*spacing(real(particles(2)%uy)) ) &
        .ge. &
        (sidebyside(2)%uy) &
            .and. &
     (particles(2)%uy &
      -2*spacing(real(particles(2)%uy)) ) &
      .le. &
       (sidebyside(2)%uy) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:106]"
      print *, "  ", "sidebyside(2)%uy (", &
 sidebyside(2)%uy, &
  ") is not", &
 particles(2)%uy,&
 "within", &
  2*spacing(real(particles(2)%uy))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%uz &
        +2*spacing(real(particles(2)%uz)) ) &
        .ge. &
        (sidebyside(2)%uz) &
            .and. &
     (particles(2)%uz &
      -2*spacing(real(particles(2)%uz)) ) &
      .le. &
       (sidebyside(2)%uz) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:107]"
      print *, "  ", "sidebyside(2)%uz (", &
 sidebyside(2)%uz, &
  ") is not", &
 particles(2)%uz,&
 "within", &
  2*spacing(real(particles(2)%uz))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(particles(1)%rod)) then
      print *, " *Assert_True failed* in test exchange &
              &[pt.fun:108]"
      print *, "  ", "particles(1)%rod is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(particles(2)%rod)) then
      print *, " *Assert_True failed* in test exchange &
              &[pt.fun:109]"
      print *, "  ", "particles(2)%rod is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

    !! check that energy is preserved in the transfer
    E_ss = potential((/sidebyside(1)%ux, sidebyside(1)%uy, sidebyside(1)%uz &
      /), (/ sidebyside(2)%ux, sidebyside(2)%uy, sidebyside(2)%uz /), &
      (/sidebyside(1)%x - sidebyside(2)%x, sidebyside(1)%y - sidebyside(2)%y,&
       sidebyside(1)%z - sidebyside(2)%z/))
    E_particles = potential((/particles(1)%ux, particles(1)%uy, &
      particles(1)%uz /), (/ particles(2)%ux, particles(2)%uy, &
      particles(2)%uz /), (/particles(1)%x - particles(2)%x, particles(1)%y &
      - particles(2)%y, particles(1)%z - particles(2)%z/))
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (E_ss &
        +2*spacing(real(E_ss)) ) &
        .ge. &
        (E_particles) &
            .and. &
     (E_ss &
      -2*spacing(real(E_ss)) ) &
      .le. &
       (E_particles) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:120]"
      print *, "  ", "E_particles (", &
 E_particles, &
  ") is not", &
 E_ss,&
 "within", &
  2*spacing(real(E_ss))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    
    !! Check that beta has been exchanged although it is not needed to be 
    !! exchanged.
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (beta_1 &
        +2*spacing(real(beta_1)) ) &
        .ge. &
        (beta) &
            .and. &
     (beta_1 &
      -2*spacing(real(beta_1)) ) &
      .le. &
       (beta) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:124]"
      print *, "  ", "beta (", &
 beta, &
  ") is not", &
 beta_1,&
 "within", &
  2*spacing(real(beta_1))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

    !! Check that energy has been exchanged
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (E_particles &
        +2*spacing(real(E_particles)) ) &
        .ge. &
        (energy) &
            .and. &
     (E_particles &
      -2*spacing(real(E_particles)) ) &
      .le. &
       (energy) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:127]"
      print *, "  ", "energy (", &
 energy, &
  ") is not", &
 E_particles,&
 "within", &
  2*spacing(real(E_particles))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

    !! check that random numbers have been exchanged
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (rand_1 &
        +2*spacing(real(rand_1)) ) &
        .ge. &
        (rand) &
            .and. &
     (rand_1 &
      -2*spacing(real(rand_1)) ) &
      .le. &
       (rand) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:130]"
      print *, "  ", "rand (", &
 rand, &
  ") is not", &
 rand_1,&
 "within", &
  2*spacing(real(rand_1))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    call mpi_barrier(MPI_COMM_WORLD, rc)    !! Makes printing sequential
  else 
    call mpi_barrier(MPI_COMM_WORLD, rc) !! Makes printing sequential

    !! check that box is now small 
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (small_box_side &
        +2*spacing(real(small_box_side)) ) &
        .ge. &
        (box_x(a_box)) &
            .and. &
     (small_box_side &
      -2*spacing(real(small_box_side)) ) &
      .le. &
       (box_x(a_box)) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:136]"
      print *, "  ", "box_x(a_box) (", &
 box_x(a_box), &
  ") is not", &
 small_box_side,&
 "within", &
  2*spacing(real(small_box_side))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (small_box_side &
        +2*spacing(real(small_box_side)) ) &
        .ge. &
        (box_y(a_box)) &
            .and. &
     (small_box_side &
      -2*spacing(real(small_box_side)) ) &
      .le. &
       (box_y(a_box)) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:137]"
      print *, "  ", "box_y(a_box) (", &
 box_y(a_box), &
  ") is not", &
 small_box_side,&
 "within", &
  2*spacing(real(small_box_side))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (small_box_side &
        +2*spacing(real(small_box_side)) ) &
        .ge. &
        (box_z(a_box)) &
            .and. &
     (small_box_side &
      -2*spacing(real(small_box_side)) ) &
      .le. &
       (box_z(a_box)) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:138]"
      print *, "  ", "box_z(a_box) (", &
 box_z(a_box), &
  ") is not", &
 small_box_side,&
 "within", &
  2*spacing(real(small_box_side))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    !! We don't test periodicity, cause it is not sent over to the other 
    !! process

    !! check that configuration is now a t-configuration
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%x &
        +2*spacing(real(particles(1)%x)) ) &
        .ge. &
        (tconf(1)%x) &
            .and. &
     (particles(1)%x &
      -2*spacing(real(particles(1)%x)) ) &
      .le. &
       (tconf(1)%x) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:143]"
      print *, "  ", "tconf(1)%x (", &
 tconf(1)%x, &
  ") is not", &
 particles(1)%x,&
 "within", &
  2*spacing(real(particles(1)%x))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%y &
        +2*spacing(real(particles(1)%y)) ) &
        .ge. &
        (tconf(1)%y) &
            .and. &
     (particles(1)%y &
      -2*spacing(real(particles(1)%y)) ) &
      .le. &
       (tconf(1)%y) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:144]"
      print *, "  ", "tconf(1)%y (", &
 tconf(1)%y, &
  ") is not", &
 particles(1)%y,&
 "within", &
  2*spacing(real(particles(1)%y))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%z &
        +2*spacing(real(particles(1)%z)) ) &
        .ge. &
        (tconf(1)%z) &
            .and. &
     (particles(1)%z &
      -2*spacing(real(particles(1)%z)) ) &
      .le. &
       (tconf(1)%z) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:145]"
      print *, "  ", "tconf(1)%z (", &
 tconf(1)%z, &
  ") is not", &
 particles(1)%z,&
 "within", &
  2*spacing(real(particles(1)%z))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%ux &
        +2*spacing(real(particles(1)%ux)) ) &
        .ge. &
        (tconf(1)%ux) &
            .and. &
     (particles(1)%ux &
      -2*spacing(real(particles(1)%ux)) ) &
      .le. &
       (tconf(1)%ux) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:146]"
      print *, "  ", "tconf(1)%ux (", &
 tconf(1)%ux, &
  ") is not", &
 particles(1)%ux,&
 "within", &
  2*spacing(real(particles(1)%ux))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%uy &
        +2*spacing(real(particles(1)%uy)) ) &
        .ge. &
        (tconf(1)%uy) &
            .and. &
     (particles(1)%uy &
      -2*spacing(real(particles(1)%uy)) ) &
      .le. &
       (tconf(1)%uy) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:147]"
      print *, "  ", "tconf(1)%uy (", &
 tconf(1)%uy, &
  ") is not", &
 particles(1)%uy,&
 "within", &
  2*spacing(real(particles(1)%uy))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(1)%uz &
        +2*spacing(real(particles(1)%uz)) ) &
        .ge. &
        (tconf(1)%uz) &
            .and. &
     (particles(1)%uz &
      -2*spacing(real(particles(1)%uz)) ) &
      .le. &
       (tconf(1)%uz) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:148]"
      print *, "  ", "tconf(1)%uz (", &
 tconf(1)%uz, &
  ") is not", &
 particles(1)%uz,&
 "within", &
  2*spacing(real(particles(1)%uz))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(particles(1)%rod)) then
      print *, " *Assert_True failed* in test exchange &
              &[pt.fun:149]"
      print *, "  ", "particles(1)%rod is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%x &
        +2*spacing(real(particles(2)%x)) ) &
        .ge. &
        (tconf(2)%x) &
            .and. &
     (particles(2)%x &
      -2*spacing(real(particles(2)%x)) ) &
      .le. &
       (tconf(2)%x) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:150]"
      print *, "  ", "tconf(2)%x (", &
 tconf(2)%x, &
  ") is not", &
 particles(2)%x,&
 "within", &
  2*spacing(real(particles(2)%x))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%y &
        +2*spacing(real(particles(2)%y)) ) &
        .ge. &
        (tconf(2)%y) &
            .and. &
     (particles(2)%y &
      -2*spacing(real(particles(2)%y)) ) &
      .le. &
       (tconf(2)%y) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:151]"
      print *, "  ", "tconf(2)%y (", &
 tconf(2)%y, &
  ") is not", &
 particles(2)%y,&
 "within", &
  2*spacing(real(particles(2)%y))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%z &
        +2*spacing(real(particles(2)%z)) ) &
        .ge. &
        (tconf(2)%z) &
            .and. &
     (particles(2)%z &
      -2*spacing(real(particles(2)%z)) ) &
      .le. &
       (tconf(2)%z) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:152]"
      print *, "  ", "tconf(2)%z (", &
 tconf(2)%z, &
  ") is not", &
 particles(2)%z,&
 "within", &
  2*spacing(real(particles(2)%z))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%ux &
        +2*spacing(real(particles(2)%ux)) ) &
        .ge. &
        (tconf(2)%ux) &
            .and. &
     (particles(2)%ux &
      -2*spacing(real(particles(2)%ux)) ) &
      .le. &
       (tconf(2)%ux) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:153]"
      print *, "  ", "tconf(2)%ux (", &
 tconf(2)%ux, &
  ") is not", &
 particles(2)%ux,&
 "within", &
  2*spacing(real(particles(2)%ux))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%uy &
        +2*spacing(real(particles(2)%uy)) ) &
        .ge. &
        (tconf(2)%uy) &
            .and. &
     (particles(2)%uy &
      -2*spacing(real(particles(2)%uy)) ) &
      .le. &
       (tconf(2)%uy) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:154]"
      print *, "  ", "tconf(2)%uy (", &
 tconf(2)%uy, &
  ") is not", &
 particles(2)%uy,&
 "within", &
  2*spacing(real(particles(2)%uy))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (particles(2)%uz &
        +2*spacing(real(particles(2)%uz)) ) &
        .ge. &
        (tconf(2)%uz) &
            .and. &
     (particles(2)%uz &
      -2*spacing(real(particles(2)%uz)) ) &
      .le. &
       (tconf(2)%uz) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:155]"
      print *, "  ", "tconf(2)%uz (", &
 tconf(2)%uz, &
  ") is not", &
 particles(2)%uz,&
 "within", &
  2*spacing(real(particles(2)%uz))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(particles(2)%rod)) then
      print *, " *Assert_True failed* in test exchange &
              &[pt.fun:156]"
      print *, "  ", "particles(2)%rod is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

    !! check that energy is preserved in the exchange
    E_t = potential((/tconf(1)%ux, tconf(1)%uy, tconf(1)%uz &
      /), (/ tconf(2)%ux, tconf(2)%uy, tconf(2)%uz /), &
      (/tconf(1)%x - tconf(2)%x, tconf(1)%y - tconf(2)%y,&
       tconf(1)%z - tconf(2)%z/))
    E_particles = potential((/particles(1)%ux, particles(1)%uy, &
      particles(1)%uz /), (/ particles(2)%ux, particles(2)%uy, &
      particles(2)%uz /), (/particles(1)%x - particles(2)%x, particles(1)%y &
      - particles(2)%y, particles(1)%z - particles(2)%z/))
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (E_t &
        +2*spacing(real(E_t)) ) &
        .ge. &
        (E_particles) &
            .and. &
     (E_t &
      -2*spacing(real(E_t)) ) &
      .le. &
       (E_particles) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:167]"
      print *, "  ", "E_particles (", &
 E_particles, &
  ") is not", &
 E_t,&
 "within", &
  2*spacing(real(E_t))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

    !! Check that energy has been exchanged
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (E_particles &
        +2*spacing(real(E_particles)) ) &
        .ge. &
        (energy) &
            .and. &
     (E_particles &
      -2*spacing(real(E_particles)) ) &
      .le. &
       (energy) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:170]"
      print *, "  ", "energy (", &
 energy, &
  ") is not", &
 E_particles,&
 "within", &
  2*spacing(real(E_particles))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

    !! Check that beta has been exchanged although it is not needed to be 
    !! exchanged.
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (beta_0 &
        +2*spacing(real(beta_0)) ) &
        .ge. &
        (beta) &
            .and. &
     (beta_0 &
      -2*spacing(real(beta_0)) ) &
      .le. &
       (beta) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:174]"
      print *, "  ", "beta (", &
 beta, &
  ") is not", &
 beta_0,&
 "within", &
  2*spacing(real(beta_0))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

    !! Check that random numbers have been exchanged. 
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (rand_0 &
        +2*spacing(real(rand_0)) ) &
        .ge. &
        (rand) &
            .and. &
     (rand_0 &
      -2*spacing(real(rand_0)) ) &
      .le. &
       (rand) )) then
      print *, " *Assert_Real_Equal failed* in test exchange &
              &[pt.fun:177]"
      print *, "  ", "rand (", &
 rand, &
  ") is not", &
 rand_0,&
 "within", &
  2*spacing(real(rand_0))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  end if

  call pt_finalize()
  call MPI_FINALIZE(rc)

  numTests = numTests + 1

 end subroutine exchange


 subroutine Setup
  noAssertFailed = .true.
 end subroutine Setup


 subroutine Teardown
 end subroutine Teardown


 subroutine test_pt( nTests, nAsserts, nAssertsTested, nFailures )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures

  continue

  call Setup
  call exchange
  call Teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_pt

end module pt_fun
