test_suite cell

!! This test is about showing the accuracy of division.
!! The test fails with tr = 1.e-16_dp.
!!
test division
  use nrtype
  implicit none
  integer :: i
  integer, parameter :: n_pos = 8
  real(dp), dimension(3, n_pos) :: positions
  type(list) :: c_list
  type(iterator) :: c_iterator
  real(dp), dimension(3) :: r
  integer, parameter :: nc = 2
  integer :: ix 
  integer :: iy
  integer :: iz
  real(dp), parameter :: tr = 1.e-15_dp

  !! Put eight positions close around the origin.
  i = 0
  do iz = 0, nc - 1
    do iy = 0, nc - 1
      do ix = 0, nc - 1
        i = i + 1
        r = (/-tr, -tr, -tr/)
        r = r + (/real(2 * ix, dp) * tr, real(2 * iy, dp) * tr, &
        real(2 * iz, dp) * tr/)
        positions(1:3, i) = r
      end do
    end do
  end do

  !! make cell list with 8 cells 
  c_list = new_list(positions, nc, nc, nc)
  do i = 1, n_pos
    c_iterator = new_iterator(c_list, i)
    !! Compare particle indices
    !! Check that all positions go to their designated cells.
    assert_equal(i, value(c_iterator))
    call advance(c_iterator)
    !! Check that there is only one position per cell. 
    assert_true(is_done(c_iterator))
  end do

  !! make cell list of only one cell
  c_list = new_list(positions, 1, 1, 1)
  !! check that all positions go to that cell
  c_iterator = new_iterator(c_list, 1)
  i = n_pos
  do while (.not. is_done(c_iterator))
    !! Compare particle indices
    assert_equal(i, value(c_iterator))
    call advance(c_iterator)
    i = i - 1
  end do
end test

test indexing
  use nrtype
  implicit none
  integer, parameter :: nx = 1
  integer, parameter :: ny = 2
  integer, parameter :: nz = 3
  type(list) :: cl
  real(dp), dimension(3, 1), parameter :: positions = reshape((/0.1_dp, 0.1_dp, 0.1_dp/), (/3, 1/))
  type(iterator) :: it
  integer :: i
  cl = new_list(positions, nx, ny, nz)
  assert_equal(0, cell_index(cl, -1, -1, -1))
  assert_equal(0, cell_index(cl, nx, ny, nz)) 

  !! Check that too small or too large cell indices produce iterator value 0
  it = new_iterator(cl, 0)
  !! Actually should test for is_done here.
  assert_equal(0, value(it))
  it = new_iterator(cl, nx * ny * nz + 1)
  assert_equal(0, value(it))

  !! Check that the one position can be found from the right cell.
  do i = 1, 6
    it = new_iterator(cl, i)
    if (.not. is_done(it)) then
      write(*, *) 'i = ', i
      exit
    end if
  end do
  assert_equal(1, value(it))
end test

end test_suite