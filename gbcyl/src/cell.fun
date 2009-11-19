test_suite cell


test one_cell_list
  use nrtype
  use particle
  use box
  implicit none
  type(boxdat) :: simbox
  real(dp) :: min_length = 10._dp
  integer :: i
  type(particledat), dimension(5) :: particles
  integer :: n_particles = 5
  type(list) :: c_list
  type(iterator) :: c_iterator
  real(dp) :: coord

  !! Create cubic box with sidelength min_length
  call make_box(simbox, min_length)
  do i = 1, n_particles
    particles(i) = new_particle()
    coord = (real(i, dp) - 0.5_dp)/real(n_particles, dp) - 0.5_dp * min_length
    call set_position(particles(i), (/coord, coord, coord/))
  end do

  !! make cell list of only one cell
  call new_list(c_list, simbox, particles, n_particles, min_length)

  !! check that all particles go to that cell
  call new_iterator(c_iterator, c_list, 1)
  i = n_particles
  do while (.not. is_done(c_iterator))
    !! Compare particle indices
    assert_equal(i, current(c_iterator))
    call next(c_iterator)
    i = i - 1
  end do
end test



test division
  use nrtype
  use particle
  use box
  implicit none
  type(boxdat) :: simbox
  real(dp) :: min_length
  integer :: i
  integer, parameter :: n_particles = 8
  type(particledat), dimension(n_particles) :: particles
  type(list) :: c_list
  type(iterator) :: c_iterator
  real(dp), dimension(3) :: r
  real(dp), parameter :: x_side = 6._dp
  real(dp) :: y_side 
  real(dp) :: z_side
  integer :: nx
  integer :: ny
  integer :: nz
  integer :: ix 
  integer :: iy
  integer :: iz
  y_side = 1.2_dp * x_side
  z_side = 1.4_dp * x_side
  call make_box(simbox, x_side)
  call set_y(simbox, y_side)
  call set_z(simbox, z_side)
  min_length = 0.499_dp * x_side
  !! make cell list 8 cells
  nx = 2
  ny = 2
  nz = 2
  do i = 1, n_particles
    particles(i) = new_particle()
    iz = (i - 1) / (nx * ny)
    iy = (i - 1 - iz * nx * ny) / nx
    ix = i - 1 - iz * nx * ny - iy * nx
    r(1) = x_side / real(nx, dp) / 2._dp + real(ix, dp) * x_side / &
    real(nx, dp) - 0.5_dp * x_side
    r(2) = y_side / real(ny, dp) / 2._dp + real(iy, dp) * y_side / &
    real(ny, dp) - 0.5_dp * y_side
    r(3) = z_side / real(nz, dp) / 2._dp + real(iz, dp) * z_side / &
    real(nz, dp) - 0.5_dp * z_side
    call set_position(particles(i), r) 
  end do 

  call new_list(c_list, simbox, particles, n_particles, min_length)
  !! check that all particles go to their designated cells
  do i = 1, n_particles
    call new_iterator(c_iterator, c_list, i)
    !! Compare particle indices
    assert_equal(i, current(c_iterator))
    call next(c_iterator)
    assert_true(is_done(c_iterator))
  end do
end test


end test_suite 