module psi6_module
use particle, only: particledat
use nrtype, only: dp
use geometry, only: minimum_image
implicit none


contains 

! Calculates the bulk bond orientational order parameter to @p particles. 
! 
! @author Jouni Karjalainen
! @author Juho Lintuvuori
!
! @p n_particles number of particles
! @p particles the array of particles
! @p Lx system size in x-direction
! @p Ly system size in y-direction
! @p Lz system size in z-direction
! @p lvec vector perpendicular to the particle planes
! 
! @return the bulk bond orientational order parameter
!
function psi6_bulk(particles, n_particles, Lx, Ly, Lz, lvec) 
  implicit none 
  double precision :: psi6_bulk
  integer, intent(in) :: n_particles
  type(particledat), dimension(n_particles), target, intent(in) :: particles
  double precision, intent(in) :: Lx, Ly, Lz
  double precision, dimension(3), intent(in) :: lvec

  integer :: i
  integer :: nrod
  complex :: psi6_bulkr
  type(particledat), pointer :: particlei

  psi6_bulk = 0.0
  psi6_bulkr = 0.0
  do i = 1, n_particles
    if (particles(i)%rod) then
      particlei => particles(i)
      psi6_bulkr = psi6_bulkr + psi6(particlei, particles, n_particles, lvec, Lx, Ly, Lz)
    end if 
  end do
  nrod = count(particles(1:n_particles)%rod)
  psi6_bulkr = psi6_bulkr / nrod
  psi6_bulk = real(psi6_bulkr)
end function psi6_bulk



! Calculates the local bond orientational order for @p particlei with respect 
! to @p particles in the direction perpendicular to @p lvec. 
! 
! @p particlei pointer to the particle to which the function is calculated
! @p particles the array of particles with respect to calculate
! @p n_particles the number of particles with respect to calculate
! @p lvec the vector perpendicular to the layers
! @p Lx system size in x-direction
! @p Ly system size in y-direction
! @p Lz system size in z-direction
!
function psi6(particlei, particles, n_particles, lvec, Lx, Ly, Lz)
  implicit none
  real(dp) :: psi6
  integer, intent(in) :: n_particles
  type(particledat), pointer :: particlei
  type(particledat), dimension(n_particles), target, intent(in) :: particles
  double precision, intent(in) :: Lx, Ly, Lz
  double precision, dimension(3), intent(in) :: lvec
  !type(geometrytype), intent(in) :: geom

  double precision :: theta
  double precision :: r 
  double precision :: pi
  complex :: arg
  double precision :: wij, sumwij
  integer :: j
  double precision, dimension(3) :: rij, u
  complex, parameter :: ii = (0.0, 1.0)

  pi = 4.0*atan(1.0)
  sumwij = 0.0
  wij = 0.0
  psi6 = 0.0
  arg = 0.0
  do j = 1, n_particles
    if(associated(particlei, target=particles(j))) cycle 
    if(.not.(particles(j)%rod)) cycle

    rij = minimum_image((/ particlei%x, particlei%y, particlei%z /), &
            (/ particles(j)%x, particles(j)%y, particles(j)%z /), &
            Lx, Ly, Lz)        
    r = sqrt(dot_product(rij, rij))
    !! Bates and Luckhurst weighting
    ! wij = 1 if r < 1.4 
    ! wij = 0 if r > 1.8 
    ! linear interpolation if 1.4 <= r <= 1.8 
    if(r > 1.8) then
      cycle
    else if(r < 1.4) then
      wij = 1.0
    else 
      wij = -(r - 1.4) / (1.8 - 1.4) + 1.0
    end if
    sumwij = sumwij + wij
    !! Remove the part perpendicular to the plane from rij   
    rij = rij - dot_product(rij, lvec) * lvec
    !! Normalize to unit length
    rij = rij / sqrt(dot_product(rij, rij))
    !! Form arbitrary unit vector in the plane. 
    u(1) = 0.0
    u(2) = 1.0
    u(3) = 0.0
    !! Remove the part perpendicular to the plane from u.
    u = u - dot_product(u, lvec) * lvec
    !! Normalize to unit length (Not originally in Juho's routine. Why?)
    u = u / sqrt(dot_product(u, u))
    !! Calculate angle between rij and u.
    theta = acos(dot_product(rij, u))
    arg = arg + wij * exp(ii * 6.0 * theta)
  end do
  if(sumwij > 0.0) then
    psi6 = arg/sumwij
  end if
end function psi6


end module psi6_module
