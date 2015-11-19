!> Implements the basic functions needed to handle data for simple, 
!! uniaxial particles such as in the Gay-Berne potential.
module particle
  use utils, only: fmt_char_dp, acceptchange
  use class_poly_box, only: poly_box, minimage
  include 'rng.inc'
  use particle_mover, only: transmove, rotate
  use iso_fortran_env, only: error_unit, output_unit, dp => REAL64
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
   contains
     procedure :: downcast_assign => particledat_downcast_assign
     procedure :: particledat_assign
     generic :: assignment(=) => particledat_assign
     procedure :: energy => particle_potential
     procedure :: nvt_update => moveparticle_2
     procedure :: pair_energy
     procedure :: to_stdout
     procedure :: move => particledat_move
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
       import particledat, poly_box, dp
       type(poly_box), intent(in) :: simbox
       type(particledat), intent(in) :: particles(:)
       integer, intent(in) :: i
       real(dp), intent(out) :: energy
       logical, intent(out) :: overlap
     end subroutine particle_energy
     
     pure subroutine single_energy(aparticle, simbox, energy, err)
       import particledat, poly_box, dp
       type(particledat), intent(in) :: aparticle
       type(poly_box), intent(in) :: simbox
       real(dp), intent(out) :: energy
       integer, intent(out) :: err
     end subroutine single_energy
  end interface

type particlearray_wrapper
   class(particledat), allocatable :: arr(:)
   logical, allocatable :: mask(:)
 contains
   procedure :: wrapper_assign
   generic :: assignment(=) => wrapper_assign
   final :: wrapper_delete
end type particlearray_wrapper

type, abstract :: pair_interaction
 contains
   procedure(pair_potential), deferred :: pair_potential
   procedure(pair_force), deferred :: pair_force
   procedure(get_cutoff), deferred :: get_cutoff
end type pair_interaction

abstract interface
   pure subroutine pair_potential(this, particlei, particlej, rij, &
        energy, err)
     import
     class(pair_interaction), intent(in) :: this
     type(particledat), intent(in) :: particlei, particlej
     real(dp), intent(in) :: rij(3)
     real(dp), intent(out) :: energy
     integer, intent(out) :: err
   end subroutine pair_potential

   pure function pair_force(this, particlei, particlej, rij) result(f)
     import
     class(pair_interaction), intent(in) :: this
     type(particledat), intent(in) :: particlei, particlej
     real(dp), intent(in) :: rij(3)
     real(dp) :: f(3)
   end function pair_force
   
   pure function get_cutoff(this) result(res)
     import pair_interaction, dp
     class(pair_interaction), intent(in) :: this
     real(dp) :: res
   end function get_cutoff
end interface

contains

  subroutine particledat_downcast_assign(this, a_particle, err)
    class(particledat), intent(inout) :: this
    class(particledat), intent(in) :: a_particle
    integer, intent(out), optional :: err
    select type (a_particle)
    type is (particledat)
       this = a_particle
    class default
       if (present(err)) err = 3
    end select
  end subroutine particledat_downcast_assign

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

!> Writes @p aparticle to the outputunit @p writeunit.
subroutine to_stdout(this)
  class(particledat), intent(in) :: this
  character(len = 2) :: typeid
  if (this%rod) then
     typeid = 'gb'
  else 
     typeid = 'lj'
  end if
  write(output_unit, '(A,1X,6(' // fmt_char_dp() // ',1X))', advance='no') &
       typeid, this%x, this%y, this%z, this%ux, &
       this%uy, this%uz 
end subroutine to_stdout


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

elemental subroutine particledat_assign(dest, src)
  class(particledat), intent(inout) :: dest
  type(particledat), intent(in) :: src
  dest%x = src%x
  dest%y = src%y
  dest%z = src%z
  dest%ux = src%ux
  dest%uy = src%uy
  dest%uz = src%uz
  dest%rod = src%rod
end subroutine particledat_assign

!> Performs a trial move of @p particles(@p i). @p simbox is the
!! simulation box in which the @p particles reside. @p genstate is the
!! random number generator state. After the move, @p dE contains the
!! change in energy of the system and @p isaccepted == .true. if the
!! move was accepted. 
pure subroutine moveparticle_2(this, genstates, simbox, temperature, nbrs, &
     pair_ia, subr_single_energy, dE, n_trials, n_accepted)
  class(particledat), intent(inout) :: this
  !! Could be type(particlearray_wrapper) if beneficial:
  class(particlearray_wrapper), intent(in) :: nbrs(:) 
  type(poly_box), intent(in) :: simbox
  real(dp), intent(in) :: temperature
  type(rngstate), intent(inout) :: genstates(0:)
  class(pair_interaction), intent(in) :: pair_ia
  procedure(single_energy) :: subr_single_energy
  real(dp), intent(out) :: dE
  integer, intent(out) :: n_trials, n_accepted

  class(particledat), allocatable :: newparticle
  class(particledat), allocatable :: oldparticle
  integer :: err
  real(dp) :: enew
  real(dp) :: eold
  logical :: isaccepted
 
  enew = 0._dp
  eold = 0._dp
  dE = 0._dp
  isaccepted = .false.
  allocate(newparticle, source=this)
  call newparticle%move(genstates(0))
  call setposition(newparticle, minimage(simbox, position(newparticle)))
  allocate(oldparticle, source=this)
  call this%downcast_assign(newparticle)
 
  call this%energy(nbrs, pair_ia, simbox, subr_single_energy, enew, err)
  
  call this%downcast_assign(oldparticle)
  if(err == 0) then 
     call this%energy(nbrs, pair_ia, simbox, subr_single_energy, &
          enew, err)
     call acceptchange(eold, enew, temperature, genstates(0), isaccepted)
     if(isaccepted) then
        call this%downcast_assign(newparticle)
        dE = enew - eold
     end if
  end if

  if (isaccepted) then
     n_accepted = 1
  else
     n_accepted = 0
  end if
  n_trials = 1
end subroutine moveparticle_2


pure subroutine particle_potential(this, nbrs, pair_ia, simbox, &
     subr_single_energy, energy, err)
  class(particledat), intent(in) :: this
  class(particlearray_wrapper), intent(in) :: nbrs(:)
  class(pair_interaction), intent(in) :: pair_ia
  type(poly_box), intent(in) :: simbox
  procedure(single_energy) :: subr_single_energy
  real(dp), intent(out) :: energy
  integer, intent(out) :: err
  real(dp) :: e
  integer :: i
  real(dp) :: rcutoff
  real(dp), allocatable :: rijs(:, :)
  logical, allocatable :: cutoff_mask(:)
  integer :: i_nbr
  energy = 0
  err = 0
  rcutoff = pair_ia%get_cutoff()
  do i_nbr = 1, size(nbrs)
     if (err /= 0) exit
     !! Select the calculated interactions by cutoff already here
     if (allocated(rijs)) deallocate(rijs)
     allocate(rijs(3, size(nbrs(i_nbr)%arr)))
     do i = 1, size(nbrs(i_nbr)%arr)
        !! :NOTE: The chosen way is better than
        !! :NOTE: rij = minimage(simbox, position(particlei), 
        !! position(particlej)) It is significantly faster to use direct
        !! references to particle coordinates.
        rijs(:, i) = minimage(simbox, nbrs(i_nbr)%arr(i)%x-this%x,& 
             nbrs(i_nbr)%arr(i)%y-this%y, nbrs(i_nbr)%arr(i)%z-this%z)

        if (all(rijs(:, i) == 0.) .and. nbrs(i_nbr)%mask(i)) then
           err = 2
           return
        end if
     end do !! ifort does not vectorize

     if (allocated(cutoff_mask)) deallocate(cutoff_mask)
     allocate(cutoff_mask(size(nbrs(i_nbr)%arr)))
     do i = 1, size(nbrs(i_nbr)%arr)
        cutoff_mask(i) = dot_product(rijs(:, i), rijs(:, i)) < rcutoff**2
     end do
     cutoff_mask = cutoff_mask .and. nbrs(i_nbr)%mask

     do i = 1, size(nbrs(i_nbr)%arr)
        if (cutoff_mask(i)) then
           call pair_ia%pair_potential(this, nbrs(i_nbr)%arr(i), rijs(:, i), &
                e, err)
           if (err /= 0) then
              energy = huge(energy)
              return
           else
              energy = energy + e
           end if
        end if
     end do
  end do
  if (err == 0) then
     call subr_single_energy(this, simbox, e, err)
     if (err == 0) energy = energy + e
  end if
end subroutine particle_potential


pure subroutine pair_energy(this, another, pair_ia, simbox, energy, err)
  class(particledat), intent(in) :: this, another
  class(pair_interaction), intent(in) :: pair_ia
  type(poly_box), intent(in) :: simbox
  real(dp), intent(out) :: energy
  integer, intent(out) :: err
  call pair_ia%pair_potential(this, another, &
       minimage(simbox, [another%x - this%x, another%y - this%y, &
       another%z - this%z]), energy, err)
end subroutine pair_energy


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
pure subroutine particledat_move(this, genstate)    
  class(particledat), intent(inout) :: this
  type(rngstate), intent(inout) :: genstate
  real(dp) :: xn, yn, zn, uxn, uyn, uzn
  call transmove(this%x, this%y, this%z, xn, yn, zn, genstate)
  this%x = xn
  this%y = yn
  this%z = zn
  if (this%rod) then
     call rotate(this%ux, this%uy, this%uz, uxn, uyn, uzn, genstate)
     this%ux = uxn
     this%uy = uyn
     this%uz = uzn
  end if
end subroutine particledat_move

  impure elemental subroutine wrapper_delete(this)
    type(particlearray_wrapper), intent(inout) :: this
    if (allocated(this%arr)) deallocate(this%arr)
    if (allocated(this%mask)) deallocate(this%mask)
  end subroutine wrapper_delete

  impure elemental subroutine wrapper_assign(this, src)
    class(particlearray_wrapper), intent(inout) :: this
    type(particlearray_wrapper), intent(in) :: src
    if (allocated(this%arr)) deallocate(this%arr)
    allocate(this%arr(size(src%arr)), source=src%arr)
    if (allocated(this%mask)) deallocate(this%mask)
    allocate(this%mask(size(src%arr)), source=src%mask)
  end subroutine wrapper_assign
  
end module particle
