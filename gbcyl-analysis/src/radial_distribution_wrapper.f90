module radial_distribution_wrapper
use nrtype, only: dp
use particle, only: particledat
implicit none


contains 



  subroutine calculate_and_print(particles, n_particles, radius, height)
  implicit none
  external gr3d
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: n_particles
  real(dp), intent(in) :: radius
  real(dp), intent(in) :: height
    real(dp) :: Lx, Ly, Lz
    logical :: x_periodic, y_periodic, z_periodic
    integer, parameter :: maxbin = 100
    real(dp), parameter :: delr = 0.1
    real(dp), dimension(maxbin) :: histogram
    x_periodic = .false.
    y_periodic = .false.
    z_periodic = .true.
    Lx = 2.0*radius
    Ly = Lx
    Lz = height
    call gr3d(particles, n_particles, Lx, Ly, Lz,&
       x_periodic, y_periodic, z_periodic, maxbin, delr, histogram)
    write(*,*) histogram
  end subroutine calculate_and_print

end module
