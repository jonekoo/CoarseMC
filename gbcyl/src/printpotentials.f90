program printpotentials
!! A program which takes as input the coordinates of two Gay-Berne particles
!! and dimensions of a cylindrical cavity and outputs the potential energies 
!! for both. Cavity wall consists of smoothly and evenly distributed Lennard-
!! Jones (LJ) particles. The molecules have two embedded LJ sites which 
!! interact with the wall. 
!! 
use nrtype
use gayberne, only: gayberne_init, potential
use utils, only : fmt_char_dp
use ljcylinder
implicit none
real(dp), dimension(3), parameter :: ui = (/0._dp, 0._dp, 1._dp/)
real(dp), dimension(3) :: uj 
real(dp), dimension(3), parameter :: ex = (/1._dp, 0._dp, 0._dp/)
real(dp), dimension(3) :: ri = (/0._dp, 0._dp, 0._dp/)
real(dp), dimension(3) :: rj = (/0._dp, 0._dp, 0._dp/)
integer, parameter :: n = 1000
integer :: i
real(dp) :: rmax = 9._dp
real(dp) :: egbgb
real(dp) :: ewall
logical :: overlap = .false.
logical :: gboverlap
real(dp) :: rcyl = 9._dp
real(dp) :: alpha = 1._dp
real(dp) :: Kw = 8._dp
real(dp) :: sigmalb = 0.5_dp
real(dp) :: sigmaub = 1.0_dp
real(dp) :: sigma
integer :: isigma
integer :: nsigma = 1 
real(dp) :: gbshift = 0._dp
real(dp) :: rhoj 
real(dp) :: kappasigma = 4.4_dp
real(dp) :: yj = 0._dp
real(dp) :: zj = 0._dp
namelist /inputnml/ Kw, sigmaub, sigmalb, nsigma, alpha, gbshift, rcyl, rmax, yj, zj
uj = ui
call gayberne_init(kappasigma, 20._dp, 1._dp, 1._dp, 1._dp, 1._dp)
read(*, nml=inputnml)
!!write(*, nml=inputnml)
rj(2) = yj
rj(3) = zj
ri = (rcyl - gbshift) * ex
do i = 1, n
  rj = rj + rmax/real(n, dp) * ex
  rhoj = sqrt(rj(1)**2 + rj(2)**2)
  if (rhoj >= rcyl) exit
  call potential(ui, uj, rj-ri, egbgb, gboverlap)
  overlap = gboverlap
  write(*, '(2('//fmt_char_dp()//',1X))', ADVANCE='NO') rcyl-rj(1), egbgb
  if (nsigma > 1) then
    do isigma = 1, nsigma
      sigma = sigmalb + (sigmaub-sigmalb)/real(nsigma-1, dp) * real(isigma-1,dp)
      ewall = ljcylinderpotential(Kw, sigma, alpha, rhoj, rcyl) 
      write(*, '(2('//fmt_char_dp()//',1X))', ADVANCE='NO') ewall, egbgb + 2._dp * ewall
    end do
    write(*, *) ''
  else if (nsigma == 1) then 
    ewall = ljcylinderpotential(Kw, sigmalb, alpha, rhoj, rcyl) 
    write(*, '(2('//fmt_char_dp()//',1X))') ewall, egbgb + 2._dp * ewall
  end if
end do
end program
