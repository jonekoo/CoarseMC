test_suite verlet

test total_pair_e
  use particle, only: particledat
  use gayberne, only: init
  use nrtype, only: dp
  use box
  use class_poly_box
  use configurations
  use class_pair_potential
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
  call initvlist(particles, n_particles, simbox, r_list)
  call pp_init(r_cutoff)
  !! Test that a simple side by side configuration produces zero total energy
  call pair_interactions(simbox, particles, pair_energy, overlap)
  assert_real_equal(0._dp, pair_energy)
  assert_false(overlap)
  !! Check that x-configuration energy at contact distance is zero
  call make_xconf(particles, n_particles)
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(simbox, particles, pair_energy, overlap)
  assert_real_equal(0._dp, pair_energy)
  assert_false(overlap)
end test

test min_image
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
  call pair_interactions(simbox, particles, xenergy, overlap)
  assert_false(overlap)
  particles(1)%x = -(box_side - sigma_0()) / 2._dp
  particles(2)%x = (box_side - sigma_0()) / 2._dp
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(simbox, particles, pair_energy, overlap)
  assert_real_equal(0._dp, pair_energy) 
  assert_false(overlap)  

  !! Test non-zero finite potential
  separation = 2._dp * sigma_0()
  box_side = 2._dp * (cutoff + tiny(box_side))
  call make_xconf(particles, n_particles)
  particles(2)%x = separation
  simbox = new_box(box_side)
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(simbox, particles, pair_energy, overlap)
  assert_real_equal(lj_potential(separation), pair_energy)
  assert_false(overlap)
  particles(1)%x = -(box_side - separation) / 2._dp
  particles(2)%x = (box_side - separation) / 2._dp
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(simbox, particles, pair_energy, overlap)
  assert_real_equal(lj_potential(separation), pair_energy) 
  assert_false(overlap)
 
  !! Test separation larger than cutoff
  box_side = 2._dp * box_side
  simbox = new_box(box_side)
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(simbox, particles, pair_energy, overlap)
  assert_real_equal(0._dp, pair_energy)
  assert_false(overlap)  
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
  call pair_interactions(simbox, particles, pair_energy, overlap)
  assert_real_equal(0._dp, pair_energy)
  assert_true(overlap)
  particles(1)%x = -(box_side - separation) / 2._dp
  particles(2)%x = (box_side - separation) / 2._dp
  call updatelist(particles, n_particles, simbox)
  call pair_interactions(simbox, particles, pair_energy, overlap)
  assert_real_equal(0._dp, pair_energy) 
  assert_true(overlap)
end test

test auto_update
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
  call pair_interactions(simbox, particles, e_pair, overlap)
  write(*, *) 'Energy of pair interactions in region 1 is ', e_pair  
  if (overlap) then
    stop 'Should not overlap in this test!' 
  end if 
  assert_true(0.0_dp > e_pair)
  call pair_interactions(simbox, particles, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(e_pair, e_single)

  ! Move to region 3
  particles(2) = region_3
  call pair_interactions(simbox, particles, e_pair, overlap)
  if (overlap) then
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(0.0_dp, e_pair)
  call pair_interactions(simbox, particles, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(e_pair, e_single)

  ! Move back to region 1
  particles(2) = region_1
  call pair_interactions(simbox, particles, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_true(0.0_dp > e_pair)
  call pair_interactions(simbox, particles, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(e_pair, e_single)

  ! Test moving from region 2 to 
  ! region 1
  ! Start from region 2
  particles(2) = region_2
  call pair_interactions(simbox, particles, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(0.0_dp, e_pair)
  call pair_interactions(simbox, particles, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(e_pair, e_single)

  ! Move to region 1
  particles(2) = region_1
  call pair_interactions(simbox, particles, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_true(0.0_dp > e_pair)
  call pair_interactions(simbox, particles, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(e_pair, e_single)

  ! and back to region 2
  particles(2) = region_2
  call pair_interactions(simbox, particles, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(0.0_dp, e_pair)
  call pair_interactions(simbox, particles, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(e_pair, e_single)

  ! Move to region 3 
  particles(2) = region_3
  call pair_interactions(simbox, particles, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  call pair_interactions(simbox, particles, e_pair, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(0.0_dp, e_pair)
  assert_real_equal(e_pair, e_single)

  ! Swap particles
  particles(1) = region_3
  particles(2) = center
  call pair_interactions(simbox, particles, particles(2), 2, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(e_pair, e_single)
  call pair_interactions(simbox, particles, particles(1), 1, &
    e_single, overlap)
  if (overlap) then 
    stop 'Should not overlap in this test!'
  end if
  assert_real_equal(e_pair, e_single)
end test

test update_symmetry
! Check that calculation of single particle potential energy
! works in the two cases:
! a. the particle moves and the list should be updated
! b. one of particles neighbour has moved and the list should be updated
end test

end test_suite