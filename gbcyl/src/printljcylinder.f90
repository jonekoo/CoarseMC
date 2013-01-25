program printljwall
use nrtype
use ljcylinder
use utils, only : fmt_char_dp
implicit none
integer :: n = 1000
integer :: i
real(dp) :: ewall
real(dp) :: rcyl = 9._dp
real(dp) :: alpha = 1._dp
real(dp) :: Kw = 8._dp
real(dp) :: sigma = 1._dp
real(dp) :: rhoj 
namelist /inputnml/ Kw, sigma, alpha, rcyl, n
read(*, nml=inputnml)
do i = 2, n
  rhoj = rcyl - real(i-1, dp)*rcyl/real(n-1, dp)
  ewall = ljcylinderpotential(Kw, sigma, alpha, rhoj, rcyl) 
  write(*, '(2('//fmt_char_dp()//',1X))') rcyl-rhoj,  ewall
end do
end program
