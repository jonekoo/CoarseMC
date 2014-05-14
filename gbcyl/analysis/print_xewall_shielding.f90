program print_xewall_shielding
use m_shielding
use m_fileunit
use utils
use class_parameterizer
implicit none
integer, parameter :: res = 600
real(dp), parameter :: lim = 3._dp
real(dp) :: tensor(3, 3) 
integer :: i
real(dp) :: pix
type(parameterizer) :: reader
character(len=80) :: idchar
integer :: argc
real(dp) :: radius, r
character(len=80) :: arg, val

argc = command_argument_count()
if (argc < 2) then
  write(*, *) 'Too few arguments!'
  call print_help()
  stop
end if  
do i = 1, argc, 2
  call get_command_argument(i, arg)
  call get_command_argument(i + 1, val)
  if (arg == "-r") then
    read(val, *) radius
  else if (arg == "-i") then
    idchar = trim(adjustl(val))
  else
    write(*, *) 'Unknown command line argument! Valid arguments are -r'
  end if
end do

reader = new_parameterizer('inputparameters.'//trim(adjustl(idchar)), &
  logfile = "print_xe_qcoupling_log."//trim(adjustl(idchar)))

call init_shielding(reader)


write(*, '("#",' // fmt_char_dp() // ',3(1X,' // fmt_char_dp() // '))') "r", "sigma_xewall_xx", "sigma_xewall_yy",&
  "sigma_xewall_zz"
pix = lim / real(res, dp) 
do i = 0, res - 1
     r = radius - real(i, dp) * pix
     tensor = xewall_shielding_local(r, radius) 
     write(*, '(4(' // fmt_char_dp() // ',1X))') r, &
       tensor(1, 1), tensor(2, 2), tensor(3, 3)
end do

call fileunit_closeallopen()

contains 

subroutine print_help()
  write(*, *) "usage:"
  write(*, *) "print_xewall_coupling -r radius -i parameterfile_index"
end subroutine

end program
