test_suite verlet

test totalpaire
  use particle, only: particledat
  use gayberne, only: gayberne_init
  use nrtype, only: dp
  use box
  use class_poly_box
  use configurations
  use class_pair_potential
  implicit none
  type(particledat), dimension(2) :: particles
  integer :: nparticles = 2
  real(dp) :: pairenergy
  logical :: overlap
  type(poly_box) :: simbox
  real(dp) :: rlist = 6.8_dp
  real(dp) :: rcutoff = 5.5_dp
  type(verletlist) :: vl
  simbox = new_box(10._dp)
  !! initialize potential
  call gayberne_init(4.4_dp, 20._dp, 1._dp, 1._dp, 1._dp, 1._dp)
  call makesidebyside(particles, nparticles)
  !! initialize verlet
  !call initvlist(particles, nparticles, simbox, rlist)
  call verlet_init(rlist, nparticles)
  vl = new_verletlist(simbox, particles)
  call pp_init(rcutoff)
  !! Test that a simple side by side configuration produces zero total energy
  call pairinteractions(vl, simbox, particles, pairenergy, overlap)
  assert_real_equal(0._dp, pairenergy)
  assert_false(overlap)
  !! Check that x-configuration energy at contact distance is zero
  call makexconf(particles, nparticles)
  call update(vl, simbox, particles)
  call pairinteractions(vl, simbox, particles, pairenergy, overlap)
  assert_real_equal(0._dp, pairenergy)
  assert_false(overlap)
  call delete(vl)
end test

test minimage
  use particle
  use nrtype
  use box
  use class_poly_box
  use gayberne
  use lj
  use configurations
  implicit none
  type(particledat), dimension(2) :: particles
  integer :: nparticles
  real(dp) :: boxside
  type(poly_box) :: simbox
  real(dp) :: pairenergy
  real(dp) :: xenergy
  logical :: overlap
  real(dp) :: cutoff = 6.8_dp
  real(dp) :: separation
  type(verletlist) :: vl
  integer :: maxneighbours = 500
  boxside = 2._dp * (cutoff + tiny(boxside))
  call makexconf(particles, nparticles)
  simbox = new_box(boxside)
  call gayberne_init(4.4_dp, 20._dp, 1._dp, 1._dp, 1._dp, 1._dp)
  call verlet_init(cutoff, maxneighbours)
  vl = new_verletlist(simbox, particles)
!  call updatelist(particles, nparticles, simbox)
  call pairinteractions(vl, simbox, particles, xenergy, overlap)
  assert_false(overlap)
  particles(1)%x = -(boxside - getsigma0()) / 2._dp
  particles(2)%x = (boxside - getsigma0()) / 2._dp
  call update(vl, simbox, particles)
  call pairinteractions(vl, simbox, particles, pairenergy, overlap)
  assert_real_equal(0._dp, pairenergy) 
  assert_false(overlap)  

  !! Test non-zero finite potential
  separation = 2._dp * getsigma0()
  boxside = 2._dp * (cutoff + tiny(boxside))
  call makexconf(particles, nparticles)
  particles(2)%x = separation
  simbox = new_box(boxside)
  call update(vl, simbox, particles)
  call pairinteractions(vl, simbox, particles, pairenergy, overlap)
  assert_real_equal(ljpotential(separation), pairenergy)
  assert_false(overlap)
  particles(1)%x = -(boxside - separation) / 2._dp
  particles(2)%x = (boxside - separation) / 2._dp
  call update(vl, simbox, particles)
  call pairinteractions(vl, simbox, particles, pairenergy, overlap)
  assert_real_equal(ljpotential(separation), pairenergy) 
  assert_false(overlap)
 
  !! Test separation larger than cutoff
  boxside = 2._dp * boxside
  simbox = new_box(boxside)
  call update(vl, simbox, particles)
  call pairinteractions(vl, simbox, particles, pairenergy, overlap)
  assert_real_equal(0._dp, pairenergy)
  assert_false(overlap)  
  call delete(vl)
end test

test overlapping
  use particle
  use nrtype
  use box
  use class_poly_box
  use gayberne
  use lj
  use configurations
  implicit none
  type(particledat), dimension(2) :: particles
  integer :: nparticles
  real(dp) :: boxside
  type(poly_box) :: simbox
  real(dp) :: pairenergy
  logical :: overlap
  real(dp) :: cutoff = 6.8_dp
  real(dp) :: separation
  type(verletlist) :: vl
  integer :: maxneighbours = 500
  separation = 0.6_dp - 1.e-9_dp
  boxside = 2._dp * (cutoff + tiny(boxside))
  call makexconf(particles, nparticles)
  particles(2)%x = separation
  simbox = new_box(boxside)
  call verlet_init(cutoff, maxneighbours)  
  vl = new_verletlist(simbox, particles)
  call pairinteractions(vl, simbox, particles, pairenergy, overlap)
  assert_real_equal(0._dp, pairenergy)
  assert_true(overlap)
  particles(1)%x = -(boxside - separation) / 2._dp
  particles(2)%x = (boxside - separation) / 2._dp
  call update(vl, simbox, particles)
  call pairinteractions(vl, simbox, particles, pairenergy, overlap)
  assert_real_equal(0._dp, pairenergy) 
  assert_true(overlap)
  call delete(vl)
end test

!test autoupdate
! A pair of particles has three different regions to be in classified by their
! separation r.
!
! 1. r < rcutoff
! 2. rcutoff <= r <= rlist
! 3. rlist < r
!
! Potential is calculated only in case 1.
! Verlet list is updated when a particle has moved more than rupdate = 
! 0.5 * (rlist - rcutoff)
!
! To completely check the verlet list we have to move particles so that they
! a. stay inside the region 
! b. move to the adjancent region
! We can not assume that the particles can not move between regions 1 and 3 
! without stopping at region 2 because that can happen with parallel tempering. 
! We have to test both: going in to a region and going out. 
! We have to test staying in the same region 

! staying in the same region: 2 * nregions test cases, moving more or less than
! rupdate
! moving between adjacent regions: 2 * 2 * 2 test cases, moving in and out of
! region 2 to 2 adjacent regions more or less than rupdate
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
!  use particle
!  use nrtype
!  use class_poly_box
!  use box
!  implicit none
!  type(particledat), dimension(2) :: particles
!  integer :: nparticles = 2
!  type(particledat) :: center, region1, region2, region3
!  type(poly_box) :: simbox
!  logical :: overlap
!  real(dp) :: rlist = 6.8_dp
!  real(dp) :: rcutoff = 5.5_dp
!  real(dp) :: r2
!  real(dp) :: rupdate
!  real(dp) :: epair
!  real(dp) :: esingle
!  center = new_particle()
!  region1 = new_particle()
!  call setposition(region1, (/1.5_dp, 0._dp, 0._dp/))
!  region2 = new_particle()
!  r2 = (rcutoff + rlist)/2._dp
!  rupdate = (rlist - rcutoff)/2._dp
!  call setposition(region2, (/r2, 0._dp, 0._dp/))
!  region3 = new_particle()
!  call setposition(region3, position(region2) + (/rupdate, 0._dp, 0._dp/))
!  particles(1) = center
!  simbox = new_box(100._dp)

  ! Test moving from region 1 to 
  ! region 3
  ! Start from region 1
!  particles(2) = region1
!  call pairinteractions(simbox, particles, epair, overlap)
!  write(*, *) 'Energy of pair interactions in region 1 is ', epair  
!  if (overlap) then
!    stop 'Should not overlap in this test!' 
!  end if 
!  assert_true(0.0_dp > epair)
!  call pairinteractions(simbox, particles, particles(2), 2, &
!    esingle, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(epair, esingle)

  ! Move to region 3
!  particles(2) = region3
!  call pairinteractions(simbox, particles, epair, overlap)
!  if (overlap) then
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(0.0_dp, epair)
!  call pairinteractions(simbox, particles, particles(2), 2, &
!    esingle, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(epair, esingle)

  ! Move back to region 1
!  particles(2) = region1
!  call pairinteractions(simbox, particles, epair, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_true(0.0_dp > epair)
!  call pairinteractions(simbox, particles, particles(2), 2, &
!    esingle, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(epair, esingle)

  ! Test moving from region 2 to 
  ! region 1
  ! Start from region 2
!  particles(2) = region2
!  call pairinteractions(simbox, particles, epair, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(0.0_dp, epair)
!  call pairinteractions(simbox, particles, particles(2), 2, &
!    esingle, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(epair, esingle)

  ! Move to region 1
!  particles(2) = region1
!  call pairinteractions(simbox, particles, epair, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_true(0.0_dp > epair)
!  call pairinteractions(simbox, particles, particles(2), 2, &
!    esingle, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(epair, esingle)

  ! and back to region 2
!  particles(2) = region2
!  call pairinteractions(simbox, particles, epair, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(0.0_dp, epair)
!  call pairinteractions(simbox, particles, particles(2), 2, &
!    esingle, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(epair, esingle)

  ! Move to region 3 
!  particles(2) = region3
!  call pairinteractions(simbox, particles, particles(2), 2, &
!    esingle, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  call pairinteractions(simbox, particles, epair, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(0.0_dp, epair)
!  assert_real_equal(epair, esingle)

  ! Swap particles
!  particles(1) = region3
!  particles(2) = center
!  call pairinteractions(simbox, particles, particles(2), 2, &
!    esingle, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(epair, esingle)
!  call pairinteractions(simbox, particles, particles(1), 1, &
!    esingle, overlap)
!  if (overlap) then 
!    stop 'Should not overlap in this test!'
!  end if
!  assert_real_equal(epair, esingle)
!end test

!test updatesymmetry
! Check that calculation of single particle potential energy
! works in the two cases:
! a. the particle moves and the list should be updated
! b. one of particles neighbour has moved and the list should be updated
!end test

end test_suite