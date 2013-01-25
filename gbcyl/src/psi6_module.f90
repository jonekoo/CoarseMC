module psi6_module
use particle
use class_poly_box
use nrtype
implicit none


contains 

! Calculates the bulk bond orientational order parameter to @p particles. 
! 
! @author Jouni Karjalainen
! @author Juho Lintuvuori
!
! @p simbox the simulation box.
! @p particles the array of particles.
! @p mask the array defining the particles to which the psi6 parameter is 
! calculated with respect to all particles. 
! @p lvec vector perpendicular to the particle planes.
! 
! @return the bulk bond orientational order parameter
!
function psi6_bulk(simbox, particles, lvec) 
  implicit none 
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), dimension(3), intent(in) :: lvec
  real(dp) :: psi6_bulk
  integer :: nparticles
  integer :: i
  !type(particledat), pointer :: particlei
  integer :: nrod
  complex(dpc) :: psi6_bulkr
  complex(dpc) :: temp
  psi6_bulkr = cmplx(0._dp, 0._dp, dpc)
  nparticles = size(particles)
  do i = 1, nparticles
    if (particles(i)%rod) then
      !particlei=>particles(i)
      temp = psi6(simbox, particles, i, lvec)
      !write(*, *) 'psi6 = ', temp
      psi6_bulkr = psi6_bulkr + temp
    end if 
  end do
  !!:TODO: Figure out if the %rod reference is correct in this case.
  !!:TODO: Are references to components of derived datatypes elemental?
  nrod = count(particles(1:nparticles)%rod)  
  psi6_bulkr = psi6_bulkr / cmplx(nrod, 0._dp, dpc)
  psi6_bulk = abs(psi6_bulkr)
end function

! Calculates the bulk bond orientational order parameter to @p particles. 
! 
! @author Jouni Karjalainen
!
! @p simbox the simulation box.
! @p particles the array of particles.
! @p mask the array defining the particles to which the psi6 parameter is 
! calculated with respect to all particles. 
! @p lvec vector perpendicular to the particle planes.
! 
! @return the bulk bond orientational order parameter for the particles(i)
! which have mask(i) =.true.
!
function psi6_bulk_masked(simbox, particles, mask, lvec) 
  implicit none 
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  logical, dimension(:), intent(in) :: mask
  real(dp), dimension(3), intent(in) :: lvec
  real(dp) :: psi6_bulk_masked
  integer :: nparticles
  integer :: i
  complex(dpc) :: psi6_bulkr
  complex(dpc) :: temp
  psi6_bulkr = cmplx(0._dp, 0._dp, dpc)
  nparticles = size(particles)
  do i = 1, nparticles
    if (particles(i)%rod .and. mask(i)) then
      temp = psi6(simbox, particles, i, lvec)
      psi6_bulkr = psi6_bulkr + temp
    end if 
  end do
  !!:TODO: Figure out if the %rod reference is correct in this case.
  !!:TODO: Are references to components of derived datatypes elemental?
  !psi6_bulkr = psi6_bulkr / cmplx(count(mask), 0._dp, dpc)
  psi6_bulk_masked = abs(psi6_bulkr)/real(count(mask), dp)
end function

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
function psi6(simbox, particles, i, lvec)
  implicit none
  !type(particledat), pointer :: particlei
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  !real(dp), intent(in) :: Lx, Ly, Lz
  real(dp), dimension(3), intent(in) :: lvec
  complex(dpc) :: psi6
  integer :: nparticles
  real(dp) :: theta
  real(dp) :: r 
  real(dp) :: pi
  complex(dpc) :: arg
  real(dp) :: wij, sumwij
  integer :: j
  real(dp), dimension(3) :: rij, u
  complex(dpc), parameter :: ii = (0._dp, 1._dp)
  nparticles = size(particles)
  pi = 4._dp*atan(1._dp)
  sumwij = 0._dp
  wij = 0._dp
  psi6 = cmplx(0._dp, 0._dp, dpc)
  arg = cmplx(0._dp, 0._dp, dpc)
  do j = 1, nparticles
    !if(associated(particlei, particles(j))) cycle
    if(.not.(particles(j)%rod)) cycle
    if (j == i) cycle
    !rij(1) = particles(i)%x - particles(j)%x
    !rij(2) = particles(i)%y - particles(j)%y
    !rij(3) = particles(i)%z - particles(j)%z
    !rij = minimum_image(rij, Lx, Ly, Lz)       
    !! :TODO: Check if there's a reason for calculating distance the way it's
    !! :TODO: done above
    rij = minimage(simbox, position(particles(j)) - position(particles(i)))
    r = sqrt(dot_product(rij, rij))
    !! Bates and Luckhurst weighting
    ! wij = 1 if r < 1.4 
    ! wij = 0 if r > 1.8 
    ! linear interpolation if 1.4 <= r <= 1.8 
    if(r > 1.8_dp) then
      cycle
    else if(r < 1.4_dp) then
      wij = 1._dp
    else 
      wij = -(r - 1.4_dp) / (1.8_dp - 1.4_dp) + 1._dp
    end if
    sumwij = sumwij + wij
    !! Remove the part perpendicular to the plane from rij   
    rij = rij - dot_product(rij, lvec) * lvec
    !! Normalize to unit length
    !! Can length of rij be 0 after removing the perpendicular part?
    rij = rij / sqrt(dot_product(rij, rij))
    !! Form arbitrary unit vector in the plane. 
    u(1) = 0._dp
    u(2) = 1._dp
    u(3) = 0._dp
    !! Remove the part perpendicular to the plane from u.
    u = u - dot_product(u, lvec) * lvec
    !! Normalize to unit length (Not originally in Juho's routine. Why?)
    u = u / sqrt(dot_product(u, u))
    !! Calculate angle between rij and u.
    theta = acos(dot_product(rij, u))
    arg = arg + cmplx(wij, 0._dp, dpc) * &
    exp(ii * cmplx(6._dp * theta, 0._dp, dpc))
  end do
  if(sumwij > 0._dp) then
    psi6 = arg / cmplx(sumwij, 0._dp, dpc)
  end if
end function psi6



end module psi6_module
