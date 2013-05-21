module distribution
use class_poly_box
use particle
use nrtype
implicit none


real(dp), save :: slice_area
real(dp), save :: direction(3) !! also slice normal

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
!! Particles i are considered as the reference particles for which the rdf is 
!! calculated. Particles j are the other particles. 
!!
subroutine distribution_func(simbox, particles, mask_i, mask_j, maxbin, delr,&
  histogram, distance_func, bin_volume_func)
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: particles(:)
  logical, intent(in) :: mask_i(size(particles)), mask_j(size(particles))
  integer, intent(in) :: maxbin
  real(dp), intent(in) :: delr
  real(dp), dimension(maxbin), intent(out) :: histogram
  interface
    function distance_func(rij)
      use nrtype
      real(dp) :: distance_func
      real(dp), intent(in) :: rij(3)
    end function
  end interface
  interface 
    function bin_volume_func(r, dr)
      use nrtype
      real(dp) :: bin_volume_func
      real(dp), intent(in) :: r
      real(dp), intent(in) :: dr
    end function
  end interface

  integer :: i, j, bin
  real(dp) :: r, rij(3)
  real(dp) :: n_ideal_gas
  real(dp) :: rlow
  real(dp) :: density_j

  histogram(1:maxbin) = 0._dp
  
  do i=1, size(particles)
     if (.not. mask_i(i)) cycle
     do j=1, size(particles)
           if (.not. mask_j(j) .or. j==i) cycle
           rij = minimage(simbox, position(particles(j)) - position(particles(i)))
           r = distance_func(rij)
           bin = int(r/delr) + 1
           !! :NOTE: this last line must be serial, or "locking" 
           if(bin <= maxbin) histogram(bin) = histogram(bin) + 1._dp
     end do
  end do
  
  !! Ideal gas density for particles j
  density_j = real(count(mask_j), dp) / volume(simbox)
  !! Normalize by n_ideal_gas = bin_volume*ideal_gas_density and N
  do bin = 1, maxbin
     rlow = real(bin - 1,dp) * delr
     n_ideal_gas = bin_volume_func(rlow, delr) * density_j
     histogram(bin) = histogram(bin) / (real(count(mask_i), dp) * n_ideal_gas)
  end do
     
end subroutine

!> Calculates the projection of rij in some direction, which has to be set
!! beforehand (module variable direction). To be used with distribution_func.
!!
!! @p rij a vector in cartesian coordinates.
!!
!! @p return distance. 
!! 
function distance_1d(rij)
  real(dp), intent(in) :: rij(3)
  real(dp) :: distance_1d
  distance_1d=abs(dot_product(direction, rij))
end function

!> Calculates the distance of two points in space. To be used with 
!! distribution_func.
!! 
!! @p rij a vector in cartesian coordinates.
!! 
!! @return distance.
!!
function distance_3d(rij)
  real(dp), intent(in) :: rij(3)
  real(dp) :: distance_3d
  distance_3d = sqrt(dot_product(rij, rij))
end function

!> Returns the volume of a spherical shell.
!!
!! @p r_inner inner radius of the shell.
!! @p thickness of the shell.
!!
function spherical_shell_volume(r_inner, thickness) result(volume)
  real(dp), intent(in) :: r_inner, thickness
  real(dp) :: volume
  volume=4._dp*pi/3._dp*((thickness + r_inner)**3-r_inner**3)
end function

!> Returns the volume of a slice with @p thickness the slice_area needs to be 
!! set beforehand. The position of the slice is ignored. 
!! 
!! @p r_inner (ignored, set to anything).
!! @p thickness the thickness of the slice.
!! 
!! @return volume of the slice
!!
function slice_volume(r_inner, thickness) result(volume)
  real(dp), intent(in) :: r_inner, thickness
  real(dp) :: volume
  volume = 2._dp * slice_area * thickness
  !! The factor 2 above comes from the fact that actually two slices, one at z
  !! and the other at -z are considered. 
end function

!! Returns the volume of a cylindrical shell. 
!!
!! @p r_inner inner radius of the shell.
!! @p r_outer outer radius of the shell.
!! @p height height of the shell.
!!
function cylinder_shell_volume(r_inner, thickness, height) result(volume)
  real(dp), intent(in) :: r_inner, thickness
  real(dp), intent(in) :: height
  real(dp) :: volume
  volume=pi*((thickness + r_inner)**2-r_inner**2)*height
end function

end module
