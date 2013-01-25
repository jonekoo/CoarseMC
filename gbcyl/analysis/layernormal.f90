module layernormal
use nrtype
use nr
use utils
use particle
use class_poly_box
implicit none

contains

subroutine globalnormal(simbox, particles, cutoff, p2, direction)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(in) :: cutoff
  real(dp), intent(out) :: p2
  real(dp), dimension(3), intent(out) :: direction
  real(sp), dimension(3, 3) :: vectors
  real(sp), dimension(3) :: values
  integer :: nrot
  real(dp), dimension(size(particles), 3) :: localnormals
  real(sp), dimension(3, 3) :: orientationtensor
  integer, dimension(1) :: ml
  integer :: i

  do i = 1, size(particles)
    localnormals(i, 1:3) = localnormal(simbox, particles, i, cutoff)
  end do
  !! 2. Calculate global layer normal as average.
  !! 2.1. Construct the tensor from local layer normals
  !! 2.2. Diagonalize tensor
  !! 2.3. Global layer normal is the eigenvector corredponding to the largest
  !! eigenvalue.
  orientationtensor = real(globaltensor(localnormals), sp)
  call jacobi(orientationtensor, values, vectors, nrot)
  !write(*, *) size(orientationtensor, 1), size(orientationtensor, 2), size(e), size(values)
 
  !call tred2(orientationtensor, values, e, .false.)
  !call tqli(values, e, vectors)
  ml = maxloc(values)
  direction = real(vectors(:, ml(1)), dp)
  p2 = real(maxval(values), dp)
end subroutine

pure function globaltensor(localnormals)
  real(dp), dimension(:, :), intent(in) :: localnormals
  real(dp), dimension(3, 3) :: globaltensor
  integer :: i
  integer :: a, b
  do i = 1, size(localnormals)/3
    forall(a = 1:3, b = 1:3) 
      globaltensor(a, b) = globaltensor(a, b) + 3._dp * localnormals(i, a) * &
      localnormals(i, b)
    end forall 
    forall(a = 1:3) globaltensor(a, a) = globaltensor(a, a) - 1._dp
  end do
  globaltensor = globaltensor/real(2 * size(localnormals) / 3, dp)
end function

function localnormal(simbox, particles, i, cutoff) result(normal)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  real(dp), intent(in) :: cutoff
  real(dp), dimension(3) :: normal
  real(sp), dimension(3, 3) :: localtensor 
  real(dp), dimension(3) :: rij, rik
  real(dp), dimension(3) :: urij, urik
  integer :: nparticles
  real(sp), dimension(3, 3) :: vectors
  real(sp), dimension(3) :: values
  integer :: a, b
  integer :: j, k
  integer :: nrot
  real(dp), dimension(3) :: cp
  integer, dimension(1) :: ml
  integer :: nneighbours
  type(particledat), dimension(size(particles)) :: neighbours
  nparticles = size(particles)
  localtensor = 0._sp
  nneighbours = 0

  !! :NOTE: PBCs have to be taken into account also here!
  !! Solution 1: Make a distance function to poly_box for PBC handling? 
  !! Solution 2: Make a distance type object for PBC handling?
  !!
  !! The box dimensions are needed anyway so adding the routine and possible 
  !! distance datatype to the box module would be a natural choice.
  !! The distance datatype is needed if we want to make elemental routine calls
  !! to calculate distances. A function returning an array of reals can't be an
  !! elemental function.
  !!
  !! Ultimately we may come to the solution that why not have a particlesystem
  !! object and ask it the distance of particles i and j? However we want to 
  !! often access properties of individual particles anyway so the facade 
  !! approach may not be beneficial. 
  !! 
  !neighbours = pack(particles, distance(position(particles(i)), &
  !position(particles(:))) < cutoff .and. )
  !nneighbours = count(-"-) - 1
  !! Find out the neighbours of particle i
  nneighbours = 0
  do j = 1, nparticles
    if (cutoff <= mindistance(simbox, position(particles(i)), &
    position(particles(j)))) cycle
    if (j == i) cycle
    nneighbours = nneighbours + 1
    neighbours(nneighbours) = particles(j) 
    call setposition(neighbours(nneighbours), minimage(simbox, &
    position(neighbours(nneighbours)) - position(particles(i))))
  end do
  do j = 1, nneighbours
    rij = position(neighbours(j))
    urij = rij / sqrt(dot_product(rij, rij))
    !! :NOTE: The following loop could probably be optimized since now j and k
    !! :NOTE: go through the same indices.
    do k = 1, nneighbours
      rik = position(neighbours(k))
      if (k == j) cycle
      urik = rik / sqrt(dot_product(rik, rik))
      cp = crossproduct(urij, urik)
      !!write(*, *) 'cp = ', cp
      forall(a = 1:3, b = 1:3) 
        localtensor(a, b) = localtensor(a, b) + real(3._dp * cp(a) * cp(b), sp)
      end forall
      forall(a = 1:3) localtensor(a, a) = localtensor(a, a) - 1._sp
    end do
  end do
  localtensor = localtensor / real(2 * nneighbours * (nneighbours - 1), sp)
  !! 1.2. Diagonalize tensor
  call jacobi(localtensor, values, vectors, nrot)    
  !! 1.3. Local layer normal is the vector corresponding to the largest 
  !! eigenvalue
  !call tred2(localtensor, values, e, .false.)
  !call tqli(values, e, vectors)
  ml = maxloc(values)
  normal = real(vectors(:, ml(1)), dp)
end function

end module
