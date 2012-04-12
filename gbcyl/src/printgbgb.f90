program printgbgb
!> A program which takes as input the coordinates of two Gay-Berne particles
!! and dimensions of a cylindrical cavity and outputs the potential energies 
!! for both. Cavity wall consists of smoothly and evenly distributed Lennard-
!! Jones (LJ) particles. The molecules have two embedded LJ sites which 
!! interact with the wall. 
!! 
use nrtype
use gayberne
use ljcylinder
use utils, only : fmt_char_dp
use class_parameterizer
use class_parameter_writer
use m_fileunit
implicit none
real(dp), dimension(3), parameter :: ex = (/1._dp, 0._dp, 0._dp/)
real(dp), dimension(3), parameter :: ey = (/0._dp, 1._dp, 0._dp/)
real(dp), dimension(3), parameter :: ez = (/0._dp, 0._dp, 1._dp/)
real(dp), dimension(3) :: rj 
integer :: i
type(parameterizer) :: reader
type(parameter_writer) :: writer
integer :: writeunit
integer :: ios

!! Iteration parameters
integer, parameter :: n = 1000
real(dp) :: rmax = 9._dp

!! Gay-Berne potential parameters
!real(dp) :: kappasigma = 4.4_dp
!real(dp) :: kappaeps = 20._dp
!real(dp) :: mu = 1._dp
!real(dp) :: nu = 1._dp
!real(dp) :: epsilon = 1._dp
!real(dp) :: sigma = 1._dp

!! Results
real(dp) :: esidebyside
real(dp) :: ecross
real(dp) :: et
real(dp) :: eendtoend
logical :: overlapsidebyside
logical :: overlapcross
logical :: overlapt
logical :: overlapendtoend

!! Read parameters
reader = new_parameterizer("ptgbcyl.in", "parameterizer.log")
call gayberne_init(reader)
call delete(reader)
do i = 1, n
  rj = real(i, dp) * rmax/real(n-1, dp) * ey
  call potential(ez, ez, rj, esidebyside, overlapsidebyside)
  call potential(ez, ex, rj, ecross, overlapcross)
  call potential(ez, ey, rj, et, overlapt)
  call potential(ey, ey, rj, eendtoend, overlapendtoend)
  write(*, '(5('//fmt_char_dp()//',1X), 4(L1,1X))') rj(2), &
  esidebyside, ecross, et, eendtoend, overlapsidebyside, overlapcross, &
  overlapt, overlapendtoend 
end do

writeunit = fileunit_getfreeunit()
open(unit=writeunit, file="printgbgb.log", POSITION="APPEND", ACTION="WRITE", &
STATUS="UNKNOWN", iostat=ios)
if(ios /= 0) then
  write(*, *) "printgbgb: Opening log file failed"
else
  writer = new_parameter_writer(writeunit)
  call gb_writeparameters(writer)
  close(writeunit)
end if
end program
