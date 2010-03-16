test_suite cell_energy


test one_cell_list
  use nrtype
  use particle
  use box
  use class_poly_box
  use cell
  implicit none
  type(poly_box) :: simbox
  real(dp) :: min_length = 10._dp
  integer :: i
  type(particledat), dimension(5) :: particles
  integer :: n_particles = 5
  type(list) :: c_list
  type(iterator) :: c_iterator
  real(dp) :: coord

  !! Create cubic box with sidelength min_length
  simbox = new_box(min_length)
  do i = 1, n_particles
    particles(i) = new_particle()
    coord = (real(i, dp) - 0.5_dp)/real(n_particles, dp) - 0.5_dp * min_length
    call set_position(particles(i), (/coord, coord, coord/))
  end do

  !! make cell list of only one cell
  c_list = new_list(simbox, particles, min_length)

  !! check that all particles go to that cell
  c_iterator = new_iterator(c_list, 1)
  i = n_particles
  do while (.not. is_done(c_iterator))
    !! Compare particle indices
    assert_equal(i, value(c_iterator))
    call advance(c_iterator)
    i = i - 1
  end do
end test

!! This test is about showing the accuracy of division for a very large 
!! system.  
!!
!! The test fails already with tr = 1.e-14 or x_side = 3000. It may be 
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
  real(dp) :: min_length
  integer :: i
  integer, parameter :: n_particles = 8
  type(particledat), dimension(n_particles) :: particles
  type(list) :: c_list
  type(iterator) :: c_iterator
  real(dp), dimension(3) :: r
  real(dp), parameter :: x_side = 2000._dp
  real(dp) :: y_side 
  real(dp) :: z_side
  integer :: nx
  integer :: ny
  integer :: nz
  integer :: ix_ 
  integer :: iy_
  integer :: iz_
  real(dp) :: tr
  min_length = 0.499_dp * x_side
  y_side = x_side + 0.5_dp * min_length
  z_side = x_side + 0.75_dp * min_length
  simbox = new_box(x_side)
  call set_y(simbox, y_side)
  call set_z(simbox, z_side)
  !! make cell list with 8 cells
  nx = 2
  ny = 2
  nz = 2
  i = 0
  do iz_ = 0, 1
    do iy_ = 0, 1
      do ix_ = 0, 1
        tr = 1.e-13_dp 
        i = i + 1
        particles(i) = new_particle()
        r = (/-tr, -tr, -tr/)
        r = r + (/real(2 * ix_, dp) * tr, real(2 * iy_, dp) * tr, &
        real(2 * iz_, dp) * tr/)
        call set_position(particles(i), r)
      end do
    end do
  end do 
  c_list = new_list(simbox, particles, min_length)
  !! check that all particles go to their designated cells
  do i = 1, n_particles
    c_iterator = new_iterator(c_list, i)
    !! Compare particle indices
    assert_equal(i, value(c_iterator))
    call advance(c_iterator)
    assert_true(is_done(c_iterator))
  end do
end test

end test_suite 