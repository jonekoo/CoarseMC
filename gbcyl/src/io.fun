test_suite io

test finding_last
  character(len = 9) :: begin = 'jalle'
  character(len = 3) :: end = 'vie'
  integer, parameter :: test_unit = 98
  character(len = *), parameter :: test_file = 'finding_last-test.txt'
  integer :: ios
  character(len = 80) :: line
  open(unit = test_unit, file = test_file, action = 'READ', status = 'OLD', iostat = ios)
  if(ios == 0) then
    call find_last(test_unit, begin, end)
    read(test_unit, *) line
    assert_equal('jalle', trim(adjustl(line)))
    read(test_unit, *) line
    assert_equal('ville', trim(adjustl(line)))
    close(test_unit)
  else
    write(*, *) 'Could not open test file.'
  end if
  open(unit = test_unit, file = test_file, action = 'READ', status = 'OLD', iostat = ios)
  if(ios == 0) then
    end = 'EOF'
    call find_last(test_unit, begin, end)
    read(test_unit, *) line
    assert_equal('jalle', trim(adjustl(line)))
    read(test_unit, *) line
    assert_equal('teno', trim(adjustl(line)))
    close(test_unit)
  else
    write(*, *) 'Could not open test file.'
  end if
  open(unit = test_unit, file = test_file, action = 'READ', status = 'OLD', iostat = ios)
  if(ios == 0) then
    begin = 'meanwhile'
    end = 'EOF'
    call find_last(test_unit, begin, end)
    read(test_unit, *) line
    assert_equal('meanwhile', trim(adjustl(line)))
    read(test_unit, *) line
    assert_equal('kalle', trim(adjustl(line)))
    close(test_unit)
  else
    write(*, *) 'Could not open test file.'
  end if
end test

!! Tests finding last configuration written in old format
!!
test finding_last_old
  use class_poly_box
  use particle
  character(len = *), parameter :: begin = '$R:'
  character(len = *), parameter :: end = 'EOF'
  integer, parameter :: test_unit = 98
  character(len = *), parameter :: test_file = 'r9-n1290-cr.gb'
  integer :: ios
  type(poly_box) :: simbox
  type(particledat), dimension(:), pointer :: particles
  integer :: n_particles
  open(unit = test_unit, file = test_file, action = 'READ', status = 'OLD', iostat = ios)
  if(ios == 0) then
    call find_last(test_unit, begin, end)
    call read_configuration(test_unit, simbox, particles, n_particles)
    close(test_unit)
  else 
    write(*, *) 'Error! test finding_last_old: Could not open file ' // test_file 
  end if
end test_suite