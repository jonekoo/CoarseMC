test_suite io

test findinglast
  character(len = 9) :: begin = 'jalle'
  character(len = 3) :: end = 'vie'
  integer, parameter :: testunit = 98
  character(len = *), parameter :: testfile = 'finding_last-test.txt'
  integer :: ios
  character(len = 80) :: line
  open(unit = testunit, file = testfile, action = 'READ', status = 'OLD', iostat = ios)
  if(ios == 0) then
    call findlast(testunit, begin, end)
    read(testunit, *) line
    assert_equal('jalle', trim(adjustl(line)))
    read(testunit, *) line
    assert_equal('ville', trim(adjustl(line)))
    close(testunit)
  else
    write(*, *) 'Could not open test file.'
  end if
  open(unit = testunit, file = testfile, action = 'READ', status = 'OLD', iostat = ios)
  if(ios == 0) then
    end = 'EOF'
    call findlast(testunit, begin, end)
    read(testunit, *) line
    assert_equal('jalle', trim(adjustl(line)))
    read(testunit, *) line
    assert_equal('teno', trim(adjustl(line)))
    close(testunit)
  else
    write(*, *) 'Could not open test file.'
  end if
  open(unit = testunit, file = testfile, action = 'READ', status = 'OLD', &
  iostat = ios)
  if(ios == 0) then
    begin = 'meanwhile'
    end = 'EOF'
    call findlast(testunit, begin, end)
    read(testunit, *) line
    assert_equal('meanwhile', trim(adjustl(line)))
    read(testunit, *) line
    assert_equal('kalle', trim(adjustl(line)))
    close(testunit)
  else
    write(*, *) 'Could not open test file.'
  end if
end test

!! Tests finding last configuration written in old format
!!
test findinglastold
  use class_poly_box
  use particle
  character(len = *), parameter :: begin = '$R:'
  character(len = *), parameter :: end = 'EOF'
  integer, parameter :: testunit = 98
  character(len = *), parameter :: testfile = 'r9-n1290-cr.gb'
  integer :: ios
  type(poly_box) :: simbox
  type(particledat), dimension(:), pointer :: particles
  integer :: nparticles
  open(unit = testunit, file = testfile, action = 'READ', status = 'OLD', iostat = ios)
  if(ios == 0) then
    call findlast(testunit, begin, end)
    call readconfiguration(testunit, simbox, particles, nparticles)
    close(testunit)
  else 
    write(*, *) 'Error! test findinglastold: Could not open file ' // testfile 
  end if
end test_suite