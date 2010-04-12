test_suite cell

!! This test is about showing the accuracy of division.
!! The test fails with tr = 1.e-16_dp.
!!
test division
  use nrtype
  implicit none
  integer :: i
  integer, parameter :: npos = 8
  real(dp), dimension(3, npos) :: positions
  type(list) :: clist
  type(iterator) :: citerator
  real(dp), dimension(3) :: r
  integer, parameter :: nc = 2
  integer :: ixpos 
  integer :: iypos
  integer :: izpos
  real(dp), parameter :: tr = 1.e-15_dp

  !! Put eight positions close around the origin.
  i = 0
  do izpos = 0, nc - 1
    do iypos = 0, nc - 1
      do ixpos = 0, nc - 1
        i = i + 1
        r = (/-tr, -tr, -tr/)
        r = r + (/real(2 * ixpos, dp) * tr, real(2 * iypos, dp) * tr, &
        real(2 * izpos, dp) * tr/)
        positions(1:3, i) = r
      end do
    end do
  end do

  !! make cell list with 8 cells 
  clist = new_list(positions, nc, nc, nc)
  do i = 1, npos
    citerator = new_iterator(clist, i)
    !! Compare particle indices
    !! Check that all positions go to their designated cells.
    assert_equal(i, value(citerator))
    call advance(citerator)
    !! Check that there is only one position per cell. 
    assert_true(isdone(citerator))
    !! Assert that advancing beyond isdone does not change anything.
    call advance(citerator)
    assert_true(isdone(citerator))
  end do

  !! make cell list of only one cell
  clist = new_list(positions, 1, 1, 1)
  !! check that all positions go to that cell
  citerator = new_iterator(clist, 1)
  i = npos
  do while (.not. isdone(citerator))
    !! Compare particle indices
    assert_equal(i, value(citerator))
    call advance(citerator)
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
  real(dp), dimension(3, 1), parameter :: positions = &
  reshape((/0.1_dp, 0.1_dp, 0.1_dp/), (/3, 1/))
  type(iterator) :: it
  integer :: i
  cl = new_list(positions, nx, ny, nz)

  !! Test that too small indices give cell zero.
  assert_equal(0, cellindex(cl, -1, 0, 0))
  assert_equal(0, cellindex(cl, 0, -1, 0))
  assert_equal(0, cellindex(cl, 0, 0, -1))
  !! Test that first cell is correct
  assert_equal(1, cellindex(cl, 0, 0, 0))

  !! Test that too large indices produce cell zero
  assert_equal(0, cellindex(cl, nx, ny - 1, nz - 1)) 
  assert_equal(0, cellindex(cl, nx - 1, ny, nz - 1)) 
  assert_equal(0, cellindex(cl, nx - 1, ny - 1, nz)) 
  assert_equal(6, cellindex(cl, nx - 1, ny - 1, nz - 1))

  !! Check that too small or too large cell indices produce iterator value 0
  it = new_iterator(cl, 0)
  assert_true(isdone(it))
  assert_equal(0, value(it))
  call advance(it) 
  assert_true(isdone(it))
  it = new_iterator(cl, nx * ny * nz + 1)
  assert_equal(0, value(it))
  assert_true(isdone(it))
  
  !! Check that the one position can be found from the right cell.
  do i = 1, 6
    it = new_iterator(cl, i)
    if (.not. isdone(it)) then
      exit
    end if
  end do
  assert_equal(1, value(it))
end test

end test_suite