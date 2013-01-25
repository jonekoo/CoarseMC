module distribution
use class_poly_box
use particle
use nrtype
implicit none

contains


!! Calculates the radial distribution function for particles of the same kind.
!! To calculate between particles of different kind, use another routine. 
!!  
!! @author Juho Lintuvuori
!! @author Jouni Karjalainen
!!
!! Based on code presented in the book Computer Simulations of Liquids by M.P.
!! Allen and D.J. Tildesley
!!
!! @p simbox the simulation box where the particles are-
!! @p particles array of particles
!! @p maxbin maximum number of abscissae in the pair correlation function
!! @p delr the distance between abscissae/bins
!! @p histogram presenting the pair correlation function
!!
!!
subroutine gr1d(simbox, particles, direction, maxbin, bin_area, bin_height, histogram)
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: particles(:)
  real(dp), intent(in) :: direction(3)
  integer, intent(in) :: maxbin
  real(dp), intent(in) :: bin_area
  real(dp), intent(in) :: bin_height
  real(dp), intent(out) :: histogram(maxbin)

  integer :: i, j, bin
  real(dp) ::  r
  real(dp) :: n_ideal_gas, density, pi
  real(dp) :: rlow
  integer :: N

  pi = 4._dp*atan(1._dp);
  histogram(1:maxbin) = 0._dp
  N=size(particles)
  do i = 1, N-1
     do j = i+1, N
       r=distance_1d(simbox, position(particles(i)), position(particles(j)), direction)
       bin = int(r/bin_height) + 1
       if(bin <= maxbin) histogram(bin) = histogram(bin) + 2._dp
     end do
  end do

  density = real(N,dp)/volume(simbox)
  !! Normalize by average density
  do bin = 1, maxbin
     rlow = real(bin-1,dp)*bin_height
     n_ideal_gas = density * bin_area * bin_height
     histogram(bin) = histogram(bin)/real(N,dp)/n_ideal_gas
  end do
     
end subroutine

!! Returns the minimum image distance between positions @p ri and @p rj 
!! projected to the unit vector @p direction.
!!
!! @p simbox the simulation box in which the distances are measured
!! @p ri position vector 
!! @p rj position vector
!! @p direction the unit vector along which the distance is measured.
!!
function distance_1d(simbox, ri, rj, direction)
  type(poly_box), intent(in) :: simbox
  real(dp), intent(in) :: ri(3), rj(3), direction(3)
  real(dp) :: distance_1d
  real(dp) :: rij(3)
  rij=minimage(simbox, rj-ri)
  distance_1d=abs(dot_product(direction, rij))
end function

!! Returns the volume of a spherical shell.
!!
!! @p r_inner inner radius of the shell.
!! @p r_outer outer radius of the shell.
!!
function spherical_shell_volume(r_inner, r_outer) result(volume)
  real(dp), intent(in) :: r_inner, r_outer
  real(dp) :: volume
  volume=4._dp*pi/3._dp*(r_outer**3-r_inner**3)
end function

!! Returns the volume of a cylindrical shell. 
!!
!! @p r_inner inner radius of the shell.
!! @p r_outer outer radius of the shell.
!! @p height height of the shell.
!!
function cylinder_shell_volume(r_inner, r_outer, height) result(volume)
  real(dp), intent(in) :: r_inner, r_outer
  real(dp), intent(in) :: height
  real(dp) :: volume
  volume=pi*(r_outer**2-r_inner**2)*height
end function

end module
