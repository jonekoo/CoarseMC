test_suite cell_energy

test onecelllist
  use nrtype
  use particle
  use box
  use class_poly_box
  use cell
  implicit none
  type(poly_box) :: simbox
  real(dp) :: minlength = 10._dp
  integer :: i
  type(particledat), dimension(5) :: particles
  integer :: nparticles = 5
  type(list) :: clist
  type(nbriterator) :: citerator
  real(dp) :: coord

  !! Create cubic box with sidelength minlength
  simbox = new_box(minlength)
  do i = 1, nparticles
    particles(i) = new_particle()
    coord = (real(i, dp) - 0.5_dp)/real(nparticles, dp) - 0.5_dp * minlength
    call setposition(particles(i), (/coord, coord, coord/))
  end do

  !! make cell list of only one cell
  call cell_energy_init(minlength, .false.)
  clist = new_celllist(simbox, particles)

  !! check that all particles go to that cell
  citerator = new_nbriterator(clist, scaledposition(simbox, particles(1)))
  i = nparticles
  do while (.not. isdone(citerator))
    !! Compare particle indices
    assert_equal(i, value(citerator))
    call advance(citerator)
    i = i - 1
  end do
end test

!! This test is about showing the accuracy of division for a very large 
!! system.  
!!
!! The test fails already with tr = 1.e-14 or xside = 3000. It may be 
!! assumed though that there should be a very large system with very 
!! many cells for the routine not to produce valid cell division. Even 
!! if the division is not precise, it does not matter, since it does not 
!! matter if the particle ends up in the cell next to the one where it 
!! should be as long as the cell sizes are large enough compared to 
!! the sum of interaction potential range with maximum translation.
!!
test division
  use nrtype
  use particle
  use box
  use class_poly_box
  use cell
  implicit none
  type(poly_box) :: simbox
  real(dp) :: minlength
  integer :: i
  integer, parameter :: nparticles = 8
  type(particledat), dimension(nparticles) :: particles
  type(list) :: clist
  type(iterator) :: citerator
  real(dp), dimension(3) :: r
  real(dp), parameter :: xside = 2000._dp
  real(dp) :: yside 
  real(dp) :: zside
  integer :: nxlist
  integer :: nylist
  integer :: nzlist
  integer :: ixlist 
  integer :: iylist
  integer :: izlist
  real(dp) :: tr
  minlength = 0.499_dp * xside
  yside = xside + 0.5_dp * minlength
  zside = xside + 0.75_dp * minlength
  simbox = new_box(xside)
  call sety(simbox, yside)
  call setz(simbox, zside)
  !! make cell list with 8 cells
  nxlist = 2
  nylist = 2
  nzlist = 2
  i = 0
  do izlist = 0, 1
    do iylist = 0, 1
      do ixlist = 0, 1
        tr = 1.e-13_dp 
        i = i + 1
        particles(i) = new_particle()
        r = (/-tr, -tr, -tr/)
        r = r + (/real(2 * ixlist, dp) * tr, real(2 * iylist, dp) * tr, &
        real(2 * izlist, dp) * tr/)
        call setposition(particles(i), r)
      end do
    end do
  end do 
  call cell_energy_init(minlength, .false.)
  clist = new_celllist(simbox, particles)
  assert_equal(ncells(clist), 8)
  !! check that all particles go to their designated cells
  do i = 1, nparticles
    citerator = new_iterator(clist, i)
    !! Compare particle indices
    assert_equal(i, value(citerator))
    call advance(citerator)
    assert_true(isdone(citerator))
  end do
end test

test nbriteration
  !! Make system of three particles. Two of them are separated by a larger 
  !! distance than the potential cutoff. One is placed in the middle of these
  !! two so that it has both particles as neighbours.
  !! xcoords = (-1.0 0.0 1.0) * 0.75 * cutoff
  use particle
  use nrtype
  use cell
  use class_poly_box
  use cylinder
  integer :: i
  type(particledat), dimension(3) :: particles
  type(list) :: cl
  real(dp), parameter :: cutoff = 5.5_dp
  type(poly_box) :: simbox
  type(nbriterator) :: it
  !do i = 1, size(particles)
  !  particles(i) = new_particle()
  !  call setposition(particles(i), (/real(i - 2, dp) * 0.75_dp * cutoff, 0._dp, 0._dp/))
  !end do
  particles(1) = new_particle()
  call setposition(particles(1), (/-0.5_dp * cutoff, -0.5_dp * cutoff, 0._dp/))
  particles(2) = new_particle()
  call setposition(particles(2), (/0.5_dp * cutoff, -0.5_dp * cutoff, 0._dp/))
  particles(3) = new_particle()
  call setposition(particles(3), (/-0.5_dp * cutoff, 0.5_dp * cutoff, 0._dp/))
  simbox = new_cylinder(2._dp * cutoff, cutoff)
  call cell_energy_init(cutoff, .false.)
  cl = new_celllist(simbox, particles)
  assert_equal(4, ncells(cl))
  it = new_nbriterator(cl, scaledposition(simbox, particles(1)))
  assert_false(isdone(it))
  assert_equal(1, value(it)) 
  assert_false(isdone(it))
!  call writetostdout(cl)
  call advance(it)
  assert_equal(2, value(it))
  call advance(it)
  assert_false(isdone(it))
  assert_equal(3, value(it))
  call advance(it)
  assert_true(isdone(it))

  !! Test periodicity in z-direction
  !! Make such a box that three cells are created in the z-direction 
  !! and only one in the x,y-plane.
  simbox = new_cylinder(cutoff, 3.3_dp * cutoff)
  call setposition(particles(1), (/0._dp, 0._dp, -1.5_dp * cutoff/))
  call setposition(particles(2), (/0._dp, 0._dp, 0._dp * cutoff/))
  call setposition(particles(3), (/0._dp, 0._dp, 1.5_dp * cutoff/))
  cl = new_celllist(simbox, particles)
  assert_equal(3, ncells(cl))
  it = new_nbriterator(cl, scaledposition(simbox, particles(1)))
  assert_false(isdone(it))
  assert_equal(3, value(it))
  call advance(it)
  assert_false(isdone(it))
  assert_equal(1, value(it))
  call advance(it)
  assert_false(isdone(it))
  assert_equal(2, value(it))
  call advance(it)
  assert_true(isdone(it))
end test

end test_suite 