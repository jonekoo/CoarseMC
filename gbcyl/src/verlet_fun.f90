! verlet_fun.f90 - a unit test suite for verlet.f90
!
! funit generated this file from verlet.fun
! at Tue Sep 22 18:00:41 +0300 2009

module verlet_fun

 use verlet

 implicit none

 logical :: noAssertFailed

 public :: test_verlet

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0



 contains

 subroutine total_pair_e

  use particle, only: particledat
!  use verlet, only: initvlist, pair_interactions
  use gayberne, only: init
  use nrtype, only: dp
  use box
  use class_poly_box
  use configurations
  implicit none
  type(particledat), dimension(2) :: particles
  integer :: n_particles = 2
  real(dp) :: pair_energy
  logical :: overlap
  type(poly_box) :: simbox
  real(dp) :: r_list = 6.8_dp
  real(dp) :: r_cutoff = 5.5_dp
  simbox = new_box(10._dp)
  !! initialize potential
  call init(4.4_dp, 20._dp, 1._dp, 1._dp, 1._dp, 1._dp)
  call make_sidebyside(particles, n_particles)
  !! initialize verlet
  call initvlist(particles, n_particles, simbox, r_list, r_cutoff)
  !! Test that a simple side by side configuration produces zero total energy
  call pair_interactions(particles, n_particles, simbox, pair_energy, overlap)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0._dp &
        +2*spacing(real(0._dp)) ) &
        .ge. &
        (pair_energy) &
            .and. &
     (0._dp &
      -2*spacing(real(0._dp)) ) &
      .le. &
       (pair_energy) )) then
      print *, " *Assert_Real_Equal failed* in test total_pair_e &
              &[verlet.fun:27]"
      print *, "  ", "pair_energy (", &
 pair_energy, &
  ") is not", &
 0._dp,&
 "within", &
  2*spacing(real(0._dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (overlap) then
      print *, " *Assert_False failed* in test total_pair_e &
              &[verlet.fun:28]"
      print *, "  ", "overlap is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  !! Check that x-configuration energy at contact distance is zero
  call make_xconf(particles, n_particles)
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(particles, n_particles, simbox, pair_energy, overlap)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0._dp &
        +2*spacing(real(0._dp)) ) &
        .ge. &
        (pair_energy) &
            .and. &
     (0._dp &
      -2*spacing(real(0._dp)) ) &
      .le. &
       (pair_energy) )) then
      print *, " *Assert_Real_Equal failed* in test total_pair_e &
              &[verlet.fun:33]"
      print *, "  ", "pair_energy (", &
 pair_energy, &
  ") is not", &
 0._dp,&
 "within", &
  2*spacing(real(0._dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (overlap) then
      print *, " *Assert_False failed* in test total_pair_e &
              &[verlet.fun:34]"
      print *, "  ", "overlap is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine total_pair_e


 subroutine min_image

  use particle
  use nrtype
  use box
  use class_poly_box
  use gayberne
  use lj
  use configurations
  implicit none
  type(particledat), dimension(2) :: particles
  integer :: n_particles
  real(dp) :: box_side
  type(poly_box) :: simbox
  real(dp) :: pair_energy
  real(dp) :: xenergy
  logical :: overlap
  real(dp) :: cutoff = 6.8_dp
  real(dp) :: separation
  box_side = 2._dp * (cutoff + tiny(box_side))
  call make_xconf(particles, n_particles)
  simbox = new_box(box_side)
  call init(4.4_dp, 20._dp, 1._dp, 1._dp, 1._dp, 1._dp)
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(particles, n_particles, simbox, xenergy, overlap)
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (overlap) then
      print *, " *Assert_False failed* in test min_image &
              &[verlet.fun:61]"
      print *, "  ", "overlap is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  particles(1)%x = -(box_side - sigma_0()) / 2._dp
  particles(2)%x = (box_side - sigma_0()) / 2._dp
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(particles, n_particles, simbox, pair_energy, overlap)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0._dp &
        +2*spacing(real(0._dp)) ) &
        .ge. &
        (pair_energy) &
            .and. &
     (0._dp &
      -2*spacing(real(0._dp)) ) &
      .le. &
       (pair_energy) )) then
      print *, " *Assert_Real_Equal failed* in test min_image &
              &[verlet.fun:66]"
      print *, "  ", "pair_energy (", &
 pair_energy, &
  ") is not", &
 0._dp,&
 "within", &
  2*spacing(real(0._dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (overlap) then
      print *, " *Assert_False failed* in test min_image &
              &[verlet.fun:67]"
      print *, "  ", "overlap is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  !! Test non-zero finite potential
  separation = 2._dp * sigma_0()
  box_side = 2._dp * (cutoff + tiny(box_side))
  call make_xconf(particles, n_particles)
  particles(2)%x = separation
  simbox = new_box(box_side)
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(particles, n_particles, simbox, pair_energy, overlap)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (lj_potential(separation) &
        +2*spacing(real(lj_potential(separation))) ) &
        .ge. &
        (pair_energy) &
            .and. &
     (lj_potential(separation) &
      -2*spacing(real(lj_potential(separation))) ) &
      .le. &
       (pair_energy) )) then
      print *, " *Assert_Real_Equal failed* in test min_image &
              &[verlet.fun:77]"
      print *, "  ", "pair_energy (", &
 pair_energy, &
  ") is not", &
 lj_potential(separation),&
 "within", &
  2*spacing(real(lj_potential(separation)))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (overlap) then
      print *, " *Assert_False failed* in test min_image &
              &[verlet.fun:78]"
      print *, "  ", "overlap is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  particles(1)%x = -(box_side - separation) / 2._dp
  particles(2)%x = (box_side - separation) / 2._dp
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(particles, n_particles, simbox, pair_energy, overlap)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (lj_potential(separation) &
        +2*spacing(real(lj_potential(separation))) ) &
        .ge. &
        (pair_energy) &
            .and. &
     (lj_potential(separation) &
      -2*spacing(real(lj_potential(separation))) ) &
      .le. &
       (pair_energy) )) then
      print *, " *Assert_Real_Equal failed* in test min_image &
              &[verlet.fun:83]"
      print *, "  ", "pair_energy (", &
 pair_energy, &
  ") is not", &
 lj_potential(separation),&
 "within", &
  2*spacing(real(lj_potential(separation)))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (overlap) then
      print *, " *Assert_False failed* in test min_image &
              &[verlet.fun:84]"
      print *, "  ", "overlap is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
 
  !! Test separation larger than cutoff
  box_side = 2._dp * box_side
  simbox = new_box(box_side)
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(particles, n_particles, simbox, pair_energy, overlap)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0._dp &
        +2*spacing(real(0._dp)) ) &
        .ge. &
        (pair_energy) &
            .and. &
     (0._dp &
      -2*spacing(real(0._dp)) ) &
      .le. &
       (pair_energy) )) then
      print *, " *Assert_Real_Equal failed* in test min_image &
              &[verlet.fun:91]"
      print *, "  ", "pair_energy (", &
 pair_energy, &
  ") is not", &
 0._dp,&
 "within", &
  2*spacing(real(0._dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (overlap) then
      print *, " *Assert_False failed* in test min_image &
              &[verlet.fun:92]"
      print *, "  ", "overlap is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine min_image


 subroutine overlapping

  use particle
  use nrtype
  use box
  use class_poly_box
  use gayberne
  use lj
  use configurations
  implicit none
  type(particledat), dimension(2) :: particles
  integer :: n_particles
  real(dp) :: box_side
  type(poly_box) :: simbox
  real(dp) :: pair_energy
  logical :: overlap
  real(dp) :: cutoff = 6.8_dp
  real(dp) :: separation
  separation = 0.6_dp - 1.e-9_dp
  box_side = 2._dp * (cutoff + tiny(box_side))
  call make_xconf(particles, n_particles)
  particles(2)%x = separation
  simbox = new_box(box_side)
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(particles, n_particles, simbox, pair_energy, overlap)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0._dp &
        +2*spacing(real(0._dp)) ) &
        .ge. &
        (pair_energy) &
            .and. &
     (0._dp &
      -2*spacing(real(0._dp)) ) &
      .le. &
       (pair_energy) )) then
      print *, " *Assert_Real_Equal failed* in test overlapping &
              &[verlet.fun:119]"
      print *, "  ", "pair_energy (", &
 pair_energy, &
  ") is not", &
 0._dp,&
 "within", &
  2*spacing(real(0._dp))
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
    if (.not.(overlap)) then
      print *, " *Assert_True failed* in test overlapping &
              &[verlet.fun:120]"
      print *, "  ", "overlap is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  particles(1)%x = -(box_side - separation) / 2._dp
  particles(2)%x = (box_side - separation) / 2._dp
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(particles, n_particles, simbox, pair_energy, overlap)
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0._dp &
        +2*spacing(real(0._dp)) ) &
        .ge. &
        (pair_energy) &
            .and. &
     (0._dp &
      -2*spacing(real(0._dp)) ) &
      .le. &
       (pair_energy) )) then
      print *, " *Assert_Real_Equal failed* in test overlapping &
              &[verlet.fun:125]"
      print *, "  ", "pair_energy (", &
 pair_energy, &
  ") is not", &
 0._dp,&
 "within", &
  2*spacing(real(0._dp))
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
    if (.not.(overlap)) then
      print *, " *Assert_True failed* in test overlapping &
              &[verlet.fun:126]"
      print *, "  ", "overlap is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine overlapping


 subroutine auto_update

! A pair of particles has three different regions to be in classified by their
! separation r.
!
! 1. r < r_cutoff
! 2. r_cutoff <= r <= r_list
! 3. r_list < r
!
! Potential is calculated only in case 1.
! Verlet list is updated when a particle has moved more than r_update = 
! 0.5 * (r_list - r_cutoff)
!
! To completely check the verlet list we have to move particles so that they
! a. stay inside the region 
! b. move to the adjancent region
! We can not assume that the particles can not move between regions 1 and 3 
! without stopping at region 2 because that can happen with parallel tempering. 
! We have to test both: going in to a region and going out. 
! We have to test staying in the same region 

! staying in the same region: 2 * n_regions test cases, moving more or less than
! r_update
! moving between adjacent regions: 2 * 2 * 2 test cases, moving in and out of
! region 2 to 2 adjacent regions more or less than r_update
! moving between non-adjacent regions: 2 cases in and out between 1 and 3 
! There are 16 test cases!!!
! Let's first forget about the cases when the region doesn't change and do the
! tests for changing region. Test cases reduced to 10.
! Let's also forget about the cases when list update is not necessary. Test
! cases reduced to 6. 
!
! Create configurations for each region so that a move from configuration to
! another results in update
! 
  use particle
  use nrtype
  use class_poly_box
  use box
  implicit none
  type(particledat), dimension(2) :: particles
  integer :: n_particles = 2
  type(particledat) :: center, region_1, region_2, region_3
  type(poly_box) :: simbox
  logical :: overlap
  real(dp) :: r_list = 6.8_dp
  real(dp) :: r_cutoff = 5.5_dp
  real(dp) :: r_2
  real(dp) :: r_update
  real(dp) :: e_pair
  real(dp) :: e_single
  center = new_particle()
  region_1 = new_particle()
  call set_position(region_1, (/1.5_dp, 0._dp, 0._dp/))
  region_2 = new_particle()
  r_2 = (r_cutoff + r_list)/2._dp
  r_update = (r_list - r_cutoff)/2._dp
  call set_position(region_2, (/r_2, 0._dp, 0._dp/))
  region_3 = new_particle()
  call set_position(region_3, position(region_2) + (/r_update, 0._dp, 0._dp/))
  particles(1) = center
  simbox = new_box(100._dp)

  ! Test moving from region 1 to 
  ! region 3
  ! Start from region 1
  particles(2) = region_1
  call pair_interactions(particles, n_particles, simbox, e_pair, overlap)
  write(*, *) 'Energy of pair interactions in region 1 is ', e_pair  
  if (overlap) then
    stop 'Should not overlap in this test!' 
  end if 
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(0.0_dp > e_pair)) then
      print *, " *Assert_True failed* in test auto_update &
              &[verlet.fun:200]"
      print *, "  ", "0.0_dp > e_pair is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  call pair_interactions(particles, n_particles, simbox, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (e_pair &
        +2*spacing(real(e_pair)) ) &
        .ge. &
        (e_single) &
            .and. &
     (e_pair &
      -2*spacing(real(e_pair)) ) &
      .le. &
       (e_single) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:206]"
      print *, "  ", "e_single (", &
 e_single, &
  ") is not", &
 e_pair,&
 "within", &
  2*spacing(real(e_pair))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  ! Move to region 3
  particles(2) = region_3
  call pair_interactions(particles, n_particles, simbox, e_pair, overlap)
  if (overlap) then
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0.0_dp &
        +2*spacing(real(0.0_dp)) ) &
        .ge. &
        (e_pair) &
            .and. &
     (0.0_dp &
      -2*spacing(real(0.0_dp)) ) &
      .le. &
       (e_pair) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:214]"
      print *, "  ", "e_pair (", &
 e_pair, &
  ") is not", &
 0.0_dp,&
 "within", &
  2*spacing(real(0.0_dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  call pair_interactions(particles, n_particles, simbox, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (e_pair &
        +2*spacing(real(e_pair)) ) &
        .ge. &
        (e_single) &
            .and. &
     (e_pair &
      -2*spacing(real(e_pair)) ) &
      .le. &
       (e_single) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:220]"
      print *, "  ", "e_single (", &
 e_single, &
  ") is not", &
 e_pair,&
 "within", &
  2*spacing(real(e_pair))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  ! Move back to region 1
  particles(2) = region_1
  call pair_interactions(particles, n_particles, simbox, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(0.0_dp > e_pair)) then
      print *, " *Assert_True failed* in test auto_update &
              &[verlet.fun:228]"
      print *, "  ", "0.0_dp > e_pair is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  call pair_interactions(particles, n_particles, simbox, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (e_pair &
        +2*spacing(real(e_pair)) ) &
        .ge. &
        (e_single) &
            .and. &
     (e_pair &
      -2*spacing(real(e_pair)) ) &
      .le. &
       (e_single) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:234]"
      print *, "  ", "e_single (", &
 e_single, &
  ") is not", &
 e_pair,&
 "within", &
  2*spacing(real(e_pair))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  ! Test moving from region 2 to 
  ! region 1
  ! Start from region 2
  particles(2) = region_2
  call pair_interactions(particles, n_particles, simbox, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0.0_dp &
        +2*spacing(real(0.0_dp)) ) &
        .ge. &
        (e_pair) &
            .and. &
     (0.0_dp &
      -2*spacing(real(0.0_dp)) ) &
      .le. &
       (e_pair) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:244]"
      print *, "  ", "e_pair (", &
 e_pair, &
  ") is not", &
 0.0_dp,&
 "within", &
  2*spacing(real(0.0_dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  call pair_interactions(particles, n_particles, simbox, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (e_pair &
        +2*spacing(real(e_pair)) ) &
        .ge. &
        (e_single) &
            .and. &
     (e_pair &
      -2*spacing(real(e_pair)) ) &
      .le. &
       (e_single) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:250]"
      print *, "  ", "e_single (", &
 e_single, &
  ") is not", &
 e_pair,&
 "within", &
  2*spacing(real(e_pair))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  ! Move to region 1
  particles(2) = region_1
  call pair_interactions(particles, n_particles, simbox, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(0.0_dp > e_pair)) then
      print *, " *Assert_True failed* in test auto_update &
              &[verlet.fun:258]"
      print *, "  ", "0.0_dp > e_pair is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  call pair_interactions(particles, n_particles, simbox, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (e_pair &
        +2*spacing(real(e_pair)) ) &
        .ge. &
        (e_single) &
            .and. &
     (e_pair &
      -2*spacing(real(e_pair)) ) &
      .le. &
       (e_single) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:264]"
      print *, "  ", "e_single (", &
 e_single, &
  ") is not", &
 e_pair,&
 "within", &
  2*spacing(real(e_pair))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  ! and back to region 2
  particles(2) = region_2
  call pair_interactions(particles, n_particles, simbox, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0.0_dp &
        +2*spacing(real(0.0_dp)) ) &
        .ge. &
        (e_pair) &
            .and. &
     (0.0_dp &
      -2*spacing(real(0.0_dp)) ) &
      .le. &
       (e_pair) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:272]"
      print *, "  ", "e_pair (", &
 e_pair, &
  ") is not", &
 0.0_dp,&
 "within", &
  2*spacing(real(0.0_dp))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  call pair_interactions(particles, n_particles, simbox, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (e_pair &
        +2*spacing(real(e_pair)) ) &
        .ge. &
        (e_single) &
            .and. &
     (e_pair &
      -2*spacing(real(e_pair)) ) &
      .le. &
       (e_single) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:278]"
      print *, "  ", "e_single (", &
 e_single, &
  ") is not", &
 e_pair,&
 "within", &
  2*spacing(real(e_pair))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  ! Move to region 3 
  particles(2) = region_3
  call pair_interactions(particles, n_particles, simbox, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  call pair_interactions(particles, n_particles, simbox, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (0.0_dp &
        +2*spacing(real(0.0_dp)) ) &
        .ge. &
        (e_pair) &
            .and. &
     (0.0_dp &
      -2*spacing(real(0.0_dp)) ) &
      .le. &
       (e_pair) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:291]"
      print *, "  ", "e_pair (", &
 e_pair, &
  ") is not", &
 0.0_dp,&
 "within", &
  2*spacing(real(0.0_dp))
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
    if (.not.( (e_pair &
        +2*spacing(real(e_pair)) ) &
        .ge. &
        (e_single) &
            .and. &
     (e_pair &
      -2*spacing(real(e_pair)) ) &
      .le. &
       (e_single) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:292]"
      print *, "  ", "e_single (", &
 e_single, &
  ") is not", &
 e_pair,&
 "within", &
  2*spacing(real(e_pair))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  ! Swap particles
  particles(1) = region_3
  particles(2) = center
  call pair_interactions(particles, n_particles, simbox, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (e_pair &
        +2*spacing(real(e_pair)) ) &
        .ge. &
        (e_single) &
            .and. &
     (e_pair &
      -2*spacing(real(e_pair)) ) &
      .le. &
       (e_single) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:302]"
      print *, "  ", "e_single (", &
 e_single, &
  ") is not", &
 e_pair,&
 "within", &
  2*spacing(real(e_pair))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  call pair_interactions(particles, n_particles, simbox, particles(1), 1, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (e_pair &
        +2*spacing(real(e_pair)) ) &
        .ge. &
        (e_single) &
            .and. &
     (e_pair &
      -2*spacing(real(e_pair)) ) &
      .le. &
       (e_single) )) then
      print *, " *Assert_Real_Equal failed* in test auto_update &
              &[verlet.fun:308]"
      print *, "  ", "e_single (", &
 e_single, &
  ") is not", &
 e_pair,&
 "within", &
  2*spacing(real(e_pair))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine auto_update


 subroutine update_symmetry

! Check that calculation of single particle potential energy
! works in the two cases:
! a. the particle moves and the list should be updated
! b. one of particles neighbour has moved and the list should be updated

  numTests = numTests + 1

 end subroutine update_symmetry


 subroutine Setup
  noAssertFailed = .true.
 end subroutine Setup


 subroutine Teardown
 end subroutine Teardown


 subroutine test_verlet( nTests, nAsserts, nAssertsTested, nFailures )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures

  continue

  call Setup
  call total_pair_e
  call Teardown

  call Setup
  call min_image
  call Teardown

  call Setup
  call overlapping
  call Teardown

  call Setup
  call auto_update
  call Teardown

  call Setup
  call update_symmetry
  call Teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_verlet

end module verlet_fun
