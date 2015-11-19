module gr3dweighted
use num_kind, only: dp
use particle
use class_poly_box
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
!! @p N number of particles
!! @p particles array of particles
!! @p Lx size of simulation cell in x-direction
!! @p Ly size of simulation cell in y-direction
!! @p Lz size of simulation cell in z-direction
!! @p maxbin maximum number of abscissae in the pair correlation function
!! @p delr the distance between abscissae/bins
!! @p histogram presenting the pair correlation function
!!
!!
subroutine gr3dweighted1(particles, simbox, maxbin, delr, orientations, &
histogram, weightedhistogram)
  type(particledat), dimension(:), intent(in) :: particles
  type(poly_box), intent(in) :: simbox
  integer, intent(in) :: maxbin
  real(dp), intent(in) :: delr
  real(dp), dimension(:, :), intent(in) :: orientations
  real(dp), dimension(maxbin), intent(out) :: histogram
  real(dp), dimension(maxbin), intent(out) :: weightedhistogram

  integer :: N
  integer :: i, j, bin
  real(dp) :: r
  real(dp) :: density, pi

  N = size(particles)
  pi = 4._dp * atan(1._dp)
  weightedhistogram(1:maxbin) = 0._dp
  histogram(1:maxbin) = 0._dp

  do i = 1, N-1
    do j = i+1, N
      r = mindistance(simbox, position(particles(i)), position(particles(j)))
      bin = int(r/delr) + 1
      if(bin <= maxbin) then !! :TODO: Could this condition be removed 
        weightedhistogram(bin) = weightedhistogram(bin) + 2._dp * &
        (3._dp/2._dp*dot_product(orientations(1:3, i), orientations(1:3, j))**2 - 0.5_dp) 
        histogram(bin) = histogram(bin) + 2._dp
      end if
      !! :TODO: Check if the definition above makes any sense! Would a 
      !! :TODO: simple cos^2 suffice.
    end do
  end do

  density = real(N, dp)/volume(simbox)
  !! Normalize by n_ideal_gas = bin_volume*ideal_gas_density and N
  do bin = 1, maxbin
     !rlow = real(bin-1, dp)*delr
     !n_ideal_gas = 4._dp*pi*density/3._dp*((rlow + delr)**3-rlow**3)
     !histogram(bin) = histogram(bin)/real(N, dp)/n_ideal_gas
     if (histogram(bin) > 0._dp) weightedhistogram(bin) = weightedhistogram(bin)/histogram(bin)
     !/real(N, dp)/n_ideal_gas
  end do
  
end subroutine

end module
