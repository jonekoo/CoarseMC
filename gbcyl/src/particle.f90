!> Implements the basic functions needed to handle data for simple, 
!! uniaxial particles such as in the Gay-Berne potential.
module particle
use nrtype, only: dp
use utils
use class_poly_box, only: poly_box, minimage
include 'rng.inc'
use particle_mover, only: transmove, rotate
implicit none

!> Holds data of a uniaxial particle, e.g. Gay-Berne 
!! particle. Using the rod variable one can use the same datatype to describe
!! two different particles.
!! 
!! x, y, z are the cartesian positional coordinates of the particle.
!! ux, uy, uz are the components of the unit vector that describes
!! the 
!! orientation of the particle.
!! rod tells if the particle is a rodlike particle or e.g. spherically
!! symmetric.
!!
!! :TODO: Remove rod variable. This module should only describe GB particles.
!!
type particledat
   real(dp) :: x = 0._dp
   real(dp) :: y = 0._dp
   real(dp) :: z = 0._dp
   real(dp) :: ux = 0._dp
   real(dp) :: uy = 0._dp
   real(dp) :: uz = 1._dp
   logical :: rod = .true.
end type particledat

interface   
   !> Returns the potential energy of @p particles(@p i).
   !! 
   !! @param simbox the simulation box where the @p particles are.
   !! @param particles the particles in the system.
   !! @param i the index of the particle in @p particles.
   !! @param energy the potential energy of @p particles(@p i)
   !! @param overlap is .true. if e.g. @p particles(@p i) is too close to
     !!        some other particle.
   !!
   pure subroutine particle_energy(simbox, particles, i, &
        energy, overlap)
     import
     type(poly_box), intent(in) :: simbox
     type(particledat), intent(in) :: particles(:)
     integer, intent(in) :: i
     real(dp), intent(out) :: energy
     logical, intent(out) :: overlap
   end subroutine particle_energy
end interface

contains

!> Writes @p aparticle to the outputunit @p writeunit.
subroutine writeparticle(writeunit, aparticle)
  integer, intent(in) :: writeunit
  type(particledat), intent(in) :: aparticle
  character(len = 2) :: typeid
  if (aparticle%rod) then
     typeid = 'gb'
  else 
     typeid = 'lj'
  end if
  write(writeunit, '(A,1X,6(' // fmt_char_dp() // ',1X))', advance='no') &
       typeid, aparticle%x, aparticle%y, aparticle%z, aparticle%ux, &
       aparticle%uy, aparticle%uz 
end subroutine writeparticle

!> Reads @p aparticle from @p readunit. If the reading succeeds, 
!! @p ios==0.
subroutine readparticle(readunit, aparticle, ios)
  type(particledat), intent(out) :: aparticle
  integer, intent(in) :: readunit
  integer, intent(out) :: ios
  character(len = 2) :: temp
  read(readunit, fmt=*, iostat=ios) &
       temp, aparticle%x, aparticle%y, aparticle%z, aparticle%ux, &
       aparticle%uy, aparticle%uz 
  if (temp == 'gb') then
     aparticle%rod = .true.
  else if (temp == 'lj') then 
     aparticle%rod = .false. 
  else
     ios = 999
  end if
  if (ios /= 0) backspace readunit
end subroutine readparticle

!> Performs a trial move of @p particles(@p i). @p simbox is the
!! simulation box in which the @p particles reside. @p genstate is the
!! random number generator state. After the move, @p dE contains the
!! change in energy of the system and @p isaccepted == .true. if the
!! move was accepted. 
pure subroutine moveparticle_2(simbox, particles, i, temperature, genstate, &
     subr_particle_energy, dE, n_trials, n_accepted)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(inout) :: particles
  integer, intent(in) :: i
  real(dp), intent(in) :: temperature
  type(rngstate), intent(inout) :: genstate
  procedure(particle_energy) :: subr_particle_energy
  real(dp), intent(out) :: dE
  integer, intent(out), optional :: n_trials, n_accepted

  type(particledat) :: newparticle
  type(particledat) :: oldparticle
  logical :: overlap
  real(dp) :: enew
  real(dp) :: eold
  logical :: isaccepted
  
  enew = 0._dp
  eold = 0._dp
  dE = 0._dp
  overlap = .false.
  isaccepted = .false.
  newparticle = particles(i)
  call move(newparticle, genstate)
  call setposition(newparticle, minimage(simbox, position(newparticle)))
  oldparticle = particles(i) 
  particles(i) = newparticle
  call subr_particle_energy(simbox, particles, i, enew, overlap)
  
  particles(i) = oldparticle
  if(.not. overlap) then 
     call subr_particle_energy(simbox, particles, i, eold, overlap)
     call acceptchange(eold, enew, temperature, genstate, isaccepted)
     if(isaccepted) then
        particles(i) = newparticle
        dE = enew - eold
     end if
  end if

  if (present(n_trials)) n_trials = 1
  if (present(n_accepted)) then
     if (isaccepted) then
        n_accepted = 1
     else
        n_accepted = 0
     end if
  end if
  
end subroutine moveparticle_2

!> Returns the position of @p aparticle as a vector.
pure function position(aparticle)
  real(dp), dimension(3) :: position
  type(particledat), intent(in) :: aparticle
  position = (/aparticle%x, aparticle%y, aparticle%z/)
end function position

!> Returns the orientation of @p aparticle as a vector. @p vec must
!! be a unit vector.
pure function orientation(aparticle) 
  real(dp), dimension(3) :: orientation
  type(particledat), intent(in) :: aparticle
  orientation = (/aparticle%ux, aparticle%uy, aparticle%uz/)
end function orientation

!> Assigns the position @p vec to @p aparticle.
pure subroutine setposition(aparticle, vec)
  type(particledat), intent(inout) :: aparticle
  real(dp), dimension(3), intent(in) :: vec
  aparticle%x = vec(1)
  aparticle%y = vec(2)
  aparticle%z = vec(3)
end subroutine setposition

!> Assigns the orientation @p vec to @p aparticle.
pure subroutine setorientation(aparticle, vec)
  type(particledat), intent(inout) :: aparticle
  real(dp), dimension(3), intent(in) :: vec
  aparticle%ux = vec(1)
  aparticle%uy = vec(2)
  aparticle%uz = vec(3)
end subroutine setorientation

!> Performs a combined translation and rotation of @p particle. 
!! @p genstate is the random number generator state.
pure subroutine move(aparticle, genstate)    
  type(particledat), intent(inout) :: aparticle
  type(rngstate), intent(inout) :: genstate
  type(particledat) :: temp
  call move2(aparticle, temp, genstate)
  aparticle = temp
end subroutine move

!> Creates a new particle @p newp by applying a combined translation
!! and rotation to @p oldp. @p genstate is the random number
!! generator state.
pure subroutine move2(oldp, newp, genstate)
  type(particledat), intent(in) :: oldp
  type(particledat), intent(out) :: newp
  type(rngstate), intent(inout) :: genstate
  newp = oldp
  call transmove(oldp%x,oldp%y,oldp%z,newp%x,newp%y,newp%z, genstate)
  if (oldp%rod) then
     call rotate(oldp%ux,oldp%uy,oldp%uz,newp%ux,newp%uy,newp%uz, genstate)
  end if
end subroutine move2
  
end module particle
