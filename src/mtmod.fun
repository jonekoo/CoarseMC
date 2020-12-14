test_suite mtmod

test restart
  integer, parameter :: mtunit = 12
  character(len = *), parameter :: mtfile = 'mt-state.txt'
  integer :: ios
  real(8) :: rn1
  real(8) :: rn2
  !! initialize rng with some seed
  call sgrnd(123456)
  !! save state
  open(unit = mtunit, file = mtfile, status = 'replace', action = 'readwrite',&
  iostat = ios)
  call mtsave(mtunit, 'f')
  !! get first rn
  rn1 = grnd()
  !! load state
  rewind(mtunit)
  call mtget(mtunit, 'f')
  !! get first rn
  rn2 = grnd()
  !! compare results
  assert_equal(rn1, rn2)
end test


end test_suite