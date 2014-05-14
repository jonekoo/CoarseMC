program print_gbxe_qcoupling
use m_quadrupole_coupling
use m_fileunit
use class_parameterizer
implicit none
integer, parameter :: res = 600
real(dp) :: xx(0:res - 1, 0:res - 1) 
real(dp) :: xz(0:res - 1, 0:res - 1) 
real(dp) :: zz(0:res - 1, 0:res - 1) 
real(dp), parameter :: lim = 3._dp
! for 0 < x < 3 * sigma_0, 0 < z < 3 * sigma_0
! print chi_zz
real(dp) :: tensor(3, 3) 
character(len=80) :: format_char
integer :: xx_unit, xz_unit, zz_unit
integer :: i, j
real(dp) :: pix
type(parameterizer) :: reader
character(len=80) :: idchar
integer :: ios

read(*, *) idchar
reader = new_parameterizer('inputparameters.'//trim(adjustl(idchar)), &
  logfile = "print_xe_qcoupling_log."//trim(adjustl(idchar)))

call qcoup_init(reader)

pix = lim / real(res, dp) 
do j = 0, res - 1
  do i = 0, res - 1
     tensor = gbxe_coupling_local((i + 0.5) * pix, (j + 0.5) * pix) 
     xx(i, j) = tensor(1, 1)
     xz(i, j) = tensor(1, 3)
     zz(i, j) = tensor(3, 3)
  end do
end do

xx_unit = fileunit_getfreeunit()
open(unit = xx_unit, file = 'gbxe_chi_xx.txt', action = 'write', status = 'replace')

xz_unit = fileunit_getfreeunit()
open(unit = xz_unit, file = 'gbxe_chi_xz.txt', action = 'write', status = 'replace')

zz_unit = fileunit_getfreeunit()
open(unit = zz_unit, file = 'gbxe_chi_zz.txt', action = 'write', status = 'replace')


write(format_char, *) '(', res, '(', trim(adjustl(fmt_char_dp())), ',1X))'
do j = 0, res - 1
  write(xx_unit, format_char, iostat = ios) xx(:, j)
  if (ios /=0 ) stop 'xx'
  write(xz_unit, format_char, iostat = ios) xz(:, j)
  if (ios /=0 ) stop 'xz'
  write(zz_unit, format_char, iostat = ios) zz(:, j)
  if (ios /=0 ) stop 'zz'
end do

call fileunit_closeallopen()

end program
