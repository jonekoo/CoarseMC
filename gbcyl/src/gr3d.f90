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
subroutine gr3d(particles, N, Lx, Ly, Lz, x_periodic, y_periodic, z_periodic,& 
  maxbin, delr, histogram)
  use nrtype, only: dp
  use particle, only: particledat
  implicit none
  integer, intent(in) :: N
  type(particledat), dimension(N), intent(in) :: particles
  real(dp), intent(in) :: Lx, Ly, Lz
  logical, intent(in) :: x_periodic, y_periodic, z_periodic
  integer, intent(in) :: maxbin
  real(dp), intent(in) :: delr
  real(dp), dimension(maxbin), intent(out) :: histogram

  integer :: i, j, bin
  real(dp) :: dx, dy, dz, r
  real(dp) :: n_ideal_gas, density, pi
  real(dp) :: rlow

  pi = 4._dp*atan(1._dp);
  histogram(1:maxbin) = 0._dp

  do i = 1, N-1
     do j = i+1, N
           dx = particles(i)%x-particles(j)%x
           dy = particles(i)%y-particles(j)%y
           dz = particles(i)%z-particles(j)%z
           
           !! Minimum image convention
           if(x_periodic) dx = dx - real(nint(dx/Lx),dp)*Lx
           if(y_periodic) dy = dy - real(nint(dy/Ly),dp)*Ly
           if(z_periodic) dz = dz - real(nint(dz/Lz),dp)*Lz
           
           r = sqrt(dx*dx + dy*dy + dz*dz)
           bin = int(r/delr) + 1
           if(bin <= maxbin) histogram(bin) = histogram(bin) + 2._dp
     end do
  end do

  density = real(N,dp)/(Lx*Ly*Lz) 
  !! Normalize by n_ideal_gas = bin_volume*ideal_gas_density and N
  do bin = 1, maxbin
     rlow = real(bin-1,dp)*delr
     n_ideal_gas = 4._dp*pi*density/3._dp*((rlow + delr)**3-rlow**3)
     histogram(bin) = histogram(bin)/real(N,dp)/n_ideal_gas
  end do
     
end subroutine gr3d
