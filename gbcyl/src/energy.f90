!> This module provides functions and subroutines for computing the
!! interaction energies for the particles in the system. The module is
!! responsible for selecting which interactions are calculated. The
!! details of the calculation such as the form of the potential
!! function of the interaction of two particles are handled by other
!! modules.
module energy
use nrtype, only: dp
use particle
use particlewall
use class_simplelist
use class_poly_box
use class_parameter_writer
use class_parameterizer
use class_pair_potential
implicit none
private
 
public :: totalparticlewallenergy
public :: energy_init
public :: energy_writeparameters
public :: get_cutoff
!public :: total_by_cell
!public :: totalenergy
public :: simple_singleparticleenergy

!> True if interactions between a wall and the particles should
!! be computed.
logical, save :: iswall

!> Cutoff radius for pair interactions. 
real(dp), save :: rcutoff = 5.5_dp

!> True when the module has been correctly initialized using
!! energy_init.
logical :: is_initialized = .false.

!interface totalenergy
!   module procedure simple_totalenergy, total_by_cell
!end interface
  
interface singleparticleenergy
   module procedure simple_singleparticleenergy
end interface

contains

!> Returns the cutoff radius for the interparticle interactions.
real(dp) function get_cutoff()
  if (.not. is_initialized) stop 'Trying to access all_pairs:get_cutoff ' // &
       'before the module is initialized.'
  get_cutoff = rcutoff
end function get_cutoff

!> Initializes this module and its dependencies.
!!
!! @param reader the object responsible for reading the parameters.
!!
subroutine energy_init(reader)
  type(parameterizer), intent(in) :: reader
  call getparameter(reader, 'is_wall_on', iswall)
  if (iswall) call particlewall_init(reader)
  call getparameter(reader, 'r_cutoff', rcutoff)
  call pp_init(reader)
  is_initialized = .true.
end subroutine
 
!> Outputs the parameters of the module and its dependencies using the
!! format and unit defined by @p writer. 
subroutine energy_writeparameters(writer)
  type(parameter_writer), intent(inout) :: writer
  call writeparameter(writer, 'is_wall_on', iswall)
  if (iswall) call particlewall_writeparameters(writer)
  call writeparameter(writer, 'r_cutoff', rcutoff)
  call pp_writeparameters(writer)
end subroutine


!> Calculates the total potential energy of @p particles with
!! respect to each other and any external fields. 
!!
!! @param simbox the simulation box where the @p particles are.
!! @param particles the array of particles to which the @p energy is
!!        calculated.
!! @param energy the total energy when the routine returns.
!! @param overlap tells if e.g. two particles are too close to each
!!        other.
!! @param n_pairs optionally returns the number of pair interactions
!!        computed.
!! 
pure subroutine simple_totalenergy(simbox, particles, energy, overlap, n_pairs)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  integer, intent(out), optional :: n_pairs
  real(dp) :: energy_pairs
  real(dp) :: energy_wall
  energy = 0._dp
  energy_wall = 0._dp
  energy_pairs = 0._dp
  overlap = .false.
  call allpairinteractions(simbox, particles, energy_pairs, overlap, n_pairs)
  if (iswall .and. (.not. overlap)) then
     call totalparticlewallenergy(simbox, particles, energy_wall, overlap)
  end if
  if (.not. overlap) then
     energy = energy_pairs + energy_wall
  end if
end subroutine simple_totalenergy

          
!> Calculates the particle-wall interaction for @p particles.
!! 
!! @param simbox the simulation box where the particles and the wall
!!        are.
!! @param particles the array of particles.
!! @param energy the total particle-wall interaction at return.
!! @param overlap is true at if a particle is too close to the wall.
!! 
pure subroutine totalparticlewallenergy(simbox, particles, energy, overlap)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap 
  integer :: i
  real(dp) :: oneprtclV 
  energy = 0._dp
  overlap = .false.
  do i = 1, size(particles)
     call particlewall_potential(particles(i), simbox, oneprtclV, overlap)
     if (overlap) return
     energy = energy + oneprtclV
  end do
end subroutine totalparticlewallenergy


!> Returns the potential energy of @p particles(@p i).
!! 
!! @param simbox the simulation box where the @p particles are.
!! @param particles the particles in the system.
!! @param i the index of the particle in @p particles.
!! @param energy the potential energy of @p particles(@p i)
!! @param overlap is .true. if e.g. @p particles(@p i) is too close to
!!        some other particle.
!!
pure subroutine simple_singleparticleenergy(simbox, particles, i, energy, overlap)
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: particles(:)
  integer, intent(in) :: i
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  real(dp) :: e_particles 
  real(dp) :: e_wall 
  e_wall = 0._dp 
  e_particles = 0._dp
  energy = 0._dp
  overlap = .false.
  if (iswall) then
     call particlewall_potential(particles(i), simbox, e_wall, overlap)
  end if
  if (.not. overlap) then
     call pairinteractions(simbox, particles, i, e_particles, overlap)
  end if
  if (.not. overlap) then 
     energy = e_particles + e_wall
  end if
end subroutine simple_singleparticleenergy


!> Computes the sum of the pair interaction energies of
!! @p particles(@p i) interacting with other particles.
!! 
!! @param simbox the simulation box for the minimum image calculation.
!! @param particles the array of particles.
!! @param i the particle to which the interactions are calculated.
!! @param energy sum of all pair interactions for particles(i) and other
!!        particles.
!! @param overlap is true if some two particles overlap. Correct energy
!!        is not guaranteed when overlap = .true. 
!! @param n_pairs optionally gives the number of pairs to which the
!!        interactions were computed.
!! 
pure subroutine pairinteractions(simbox, particles, i, energy, overlap, &
     n_pairs)
  type(particledat), dimension(:), intent(in) :: particles
  type(poly_box), intent(in) :: simbox
  integer, intent(in) :: i
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  integer, intent(out), optional :: n_pairs
  integer :: nparticles
  integer :: j
  real(dp) :: epair
  logical :: cutoff_mask(size(particles))
  real(dp) :: rijs(3, size(particles))
  logical :: overlap_ij
  energy = 0._dp
  overlap = .false.
  overlap_ij = .false.
  nparticles = size(particles)
  if (present(n_pairs)) n_pairs = 0

  !! Select the calculated interactions by cutoff already here
  do j = 1, nparticles
     !! :NOTE: The chosen way is better than
     !! :NOTE: rij = minimage(simbox, position(particlei), 
     !! position(particlej)) It is significantly faster to use direct
     !! references to particle coordinates.
     rijs(:, j) = minimage(simbox, (/particles(j)%x-particles(i)%x,& 
          particles(j)%y-particles(i)%y, particles(j)%z-particles(i)%z/))
  end do !! ifort does not vectorize
  do j = 1, nparticles
     cutoff_mask(j) = dot_product(rijs(:,j), rijs(:,j)) < rcutoff**2
  end do
  !! Remove the particle i from the cutoff mask 
  cutoff_mask(i) = .false. 
  
  !! :TODO: test overlap conditions also here?
  !! :TODO: sort particles by type here as well?  
  do j = 1, nparticles
     if (cutoff_mask(j)) then
        call pair_potential(particles(i), particles(j), rijs(:,j), epair, &
             overlap)
        if (present(n_pairs)) n_pairs = n_pairs + 1
        if(overlap) then
           return
        else
           energy = energy + epair
        end if
     end if
  end do 
end subroutine pairinteractions

!> Computes all pair interactions without any neighbourlist. Used for
!! debugging.
!!
!! @param simbox the simulation box in which the particles reside.
!! @param particles the array of particles.
!! @param energy the contribution to total energy from all pair
!!        interactions.
!! @param overlap is true if e.g. two particles are too close to each
!!        other.
!! @param n_pairs optionally gives the number of pair interactions
!!        computed.
!!
pure subroutine allpairinteractions(simbox, particles, energy, overlap, &
     n_pairs)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  integer, intent(out), optional :: n_pairs
  integer :: nparticles
  integer :: i, j
  real(dp) :: epair
  real(dp) :: rij(3)
  if (present(n_pairs)) n_pairs = 0 
  energy = 0._dp
  overlap = .false.
  nparticles = size(particles)
  do i = 1, nparticles - 1
     do j = i + 1, nparticles
        rij = minimage(simbox, position(particles(j))-position(particles(i)))
        if (dot_product(rij, rij) < rcutoff**2) then
           call pair_potential(particles(i), particles(j), rij, epair, overlap)
           if (present(n_pairs)) n_pairs = n_pairs + 1
           if(overlap) then
              return
           else
              energy = energy + epair
           end if
        end if
     end do
  end do
end subroutine allpairinteractions


end module
