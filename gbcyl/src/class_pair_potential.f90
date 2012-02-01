module class_pair_potential
use nrtype
use gayberne
use particle
use class_poly_box
use class_parameterizer
use class_parameter_writer
implicit none
private

public :: pairv
public :: cutoff
public :: pp_init
public :: pp_writeparameters

real(dp), save :: rcutoff = 5.5_dp

type pairpotential
  private
  real(dp) :: cutoff
end type

interface pp_init
  module procedure pp_initwith, pp_initvalue
end interface

contains

!! Initializes the verlet neighbour list.
!! 
!! @p particles the array of particles
!! @p nparticles the number of particles
!! @p simbox the simulation cell
!! @p reader gives the parameters
!! 
subroutine pp_initwith(reader)
  type(parameterizer), intent(in) :: reader
  call getparameter(reader, 'r_cutoff', rcutoff)
  call gayberne_init(reader)
end subroutine

subroutine pp_initvalue(cutoff)
  real(dp), intent(in) :: cutoff
  rcutoff = cutoff
end subroutine

!elemental function new_pairpotential(cutoff) result(pp)
!  real(dp), intent(in) :: cutoff
!  type(pairpotential) :: pp
!  pp%cutoff = cutoff
!end function

subroutine pp_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call writeparameter(writer, 'r_cutoff', rcutoff)
  call gb_writeparameters(writer)
end subroutine

elemental function cutoff() result(r)
  real(dp) :: r
  r = rcutoff
end function

!! Calculates the interaction energy of a pair of particles. 
!! 
!! @p particlei the first particle
!! @p particlej the second particle
!! @p simbox the simulation cell
!! @p potE the interaction energy
!! @p overlap is true if the two particles overlap each other
!! 
!! :TODO: put this somewhere else. Problem is that you need to initialize 
!! the potential module to use this. 
!!
!! A poly_particle class could be used. This class could be quered for
!! the true type of the particles since it knows it anyway.
!! 
elemental subroutine pairv(particlei, particlej, simbox, potE, overlap)
  type(particledat), intent(in) :: particlei 
  type(particledat), intent(in) :: particlej
  type(poly_box), intent(in) :: simbox
  real(dp), intent(out) :: potE
  logical, intent(out) :: overlap
  real(dp), dimension(3) :: rij, ui, uj
  !potE = 0._dp
  !overlap = .false.
  !! :NOTE: The chosen way is better than
  !! :NOTE: rij = minimage(simbox, position(particlei), position(particlej))
  !! It is notably faster to use direct references to particle coordinates.
  rij = minimage(simbox, (/particlej%x-particlei%x,particlej%y-particlei%y,&
  & particlej%z-particlei%z/))
  if (dot_product(rij, rij) < rcutoff**2) then
    ui(1)=particlei%ux
    ui(2)=particlei%uy
    ui(3)=particlei%uz
    uj(1)=particlej%ux
    uj(2)=particlej%uy
    uj(3)=particlej%uz
    call potential(ui, uj, rij, potE, overlap)  
  else
    potE = 0._dp
    overlap = .false.
  end if
end subroutine

end module
