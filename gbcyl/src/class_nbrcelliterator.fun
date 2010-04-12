test_suite class_nbrcelliterator

test iterator1d
  type(nbrcelliterator1d) :: it

  !! Test for no neighbours
  it = new_nbrcelliterator1d(1, 0)
  assert_false(isdone(it))
  assert_equal(0, value(it))
  call advance(it)
  assert_true(isdone(it))

  !! Test for one neighbour only
  it = new_nbrcelliterator1d(2, 1)
  assert_false(isdone(it))
  assert_equal(value(it), 1)
  call advance(it)
  assert_false(isdone(it))
  assert_equal(value(it), 0)
  call advance(it)
  assert_true(isdone(it))
  !! Make sure that advancing beyond end does not change anything.
  call advance(it)
  assert_true(isdone(it))

  !! Test for two neighbours
  it = new_nbrcelliterator1d(3, 1)
  assert_false(isdone(it))
  assert_equal(value(it), 0)
  call advance(it)
  assert_false(isdone(it))
  assert_equal(value(it), 1)
  call advance(it)
  assert_false(isdone(it))
  assert_equal(value(it), 2)
  call advance(it)
  assert_true(isdone(it))

  !! Test periodicity
  it = new_nbrcelliterator1d(3, 0)
  assert_false(isdone(it))
  assert_equal(2, value(it))
  call advance(it)
  assert_false(isdone(it))
  assert_equal(0, value(it))
  call advance(it)
  assert_false(isdone(it))
  assert_equal(1, value(it))
  call advance(it)
  assert_true(isdone(it))

  !! Test for invalid starting points out of range
  it = new_nbrcelliterator1d(3, 3)
  assert_true(isdone(it))
  it = new_nbrcelliterator1d(3, -1)
  assert_true(isdone(it))
end test

test nbrcelliteration
  type(nbrcelliterator) :: it
  integer, parameter :: nx = 1
  integer, parameter :: ny = 2
  integer, parameter :: nz = 3
  integer :: ix, iy, iz
  ix = 0
  iy = 0
  iz = 1

  !! Test normal operation
  it = new_nbrcelliterator(nx, ny, nz, ix, iy, iz)
  assert_false(isdone(it))
  assert_equal(nbrcellit_xvalue(it), 0)
  assert_equal(nbrcellit_yvalue(it), 0)
  assert_equal(nbrcellit_zvalue(it), 0)   
  call advance(it)
  assert_false(isdone(it))
  assert_equal(nbrcellit_xvalue(it), 0)
  assert_equal(nbrcellit_yvalue(it), 1)
  assert_equal(nbrcellit_zvalue(it), 0)   
  call advance(it)
  assert_false(isdone(it))
  assert_equal(nbrcellit_xvalue(it), 0)
  assert_equal(nbrcellit_yvalue(it), 0)
  assert_equal(nbrcellit_zvalue(it), 1)   
  call advance(it)
  assert_false(isdone(it))
  assert_equal(nbrcellit_xvalue(it), 0)
  assert_equal(nbrcellit_yvalue(it), 1)
  assert_equal(nbrcellit_zvalue(it), 1)   
  call advance(it)
  assert_false(isdone(it))
  assert_equal(nbrcellit_xvalue(it), 0)
  assert_equal(nbrcellit_yvalue(it), 0)
  assert_equal(nbrcellit_zvalue(it), 2)   
  call advance(it)
  assert_false(isdone(it))
  assert_equal(nbrcellit_xvalue(it), 0)
  assert_equal(nbrcellit_yvalue(it), 1)
  assert_equal(nbrcellit_zvalue(it), 2)   
  call advance(it)
  assert_true(isdone(it))

  !! Test for out of range starting indices
  it = new_nbrcelliterator(nx, ny, nz, -1, 0, 0)
  assert_true(isdone(it))
  it = new_nbrcelliterator(nx, ny, nz, 1, 0, 0)
  assert_true(isdone(it))
  it = new_nbrcelliterator(nx, ny, nz, 0, -1, 0)
  assert_true(isdone(it))
  it = new_nbrcelliterator(nx, ny, nz, 0, 2, 0)
  assert_true(isdone(it))
  it = new_nbrcelliterator(nx, ny, nz, 0, 0, -1)
  assert_true(isdone(it))
  it = new_nbrcelliterator(nx, ny, nz, 0, 0, 3)
  assert_true(isdone(it))


  !! Test for periodicity
  it = new_nbrcelliterator(1, 1, 3, 0, 0, 0)
  assert_false(isdone(it))
  assert_equal(2, value(it))
  call advance(it)
  assert_false(isdone(it))
  assert_equal(0, value(it))
  call advance(it)
  assert_false(isdone(it))
  assert_equal(1, value(it))
  call advance(it)
  assert_true(isdone(it))
end test

end test_suite