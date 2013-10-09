module class_pair_potential
use nrtype
use gayberne
use lj, only: lj_potential, lj_init, lj_writeparameters, lj_force
use gblj, only: gblj_potential, gblj_init, gblj_force
use particle
use class_poly_box
use class_parameterizer
use class_parameter_writer
implicit none
private

public :: pairv
public :: pp_init
public :: pp_writeparameters
public :: pair_force

real(dp), save :: rcutoff = 5.5_dp

type pairpotential
  private
  real(dp) :: cutoff
end type

interface pp_init
  module procedure pp_initwith
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
  call gayberne_init(reader)
  call lj_init(reader)
end subroutine

subroutine pp_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call gb_writeparameters(writer)
  call lj_writeparameters(writer)
end subroutine

!! Calculates the interaction energy of a pair of particles. 
!! 
!! @p particlei the first particle
!! @p particlej the second particle
!! @p simbox the simulation cell
!! @p rij is the (minimum image) vector from particle i to particle j
!! @p potE the interaction energy
!! @p overlap is true if the two particles overlap each other
!! 
!! :TODO: put this somewhere else. Problem is that you need to initialize 
!! the potential module to use this. 
!!
!! A poly_particle class could be used. This class could be quered for
!! the true type of the particles since it knows it anyway.
!! 
pure subroutine pairv(particlei, particlej, rij, potE, overlap)
  type(particledat), intent(in) :: particlei 
  type(particledat), intent(in) :: particlej
  real(dp), intent(in) :: rij(3)
  real(dp), intent(out) :: potE
  logical, intent(out) :: overlap
  real(dp), dimension(3) :: ui, uj
  ui(1)=particlei%ux
  ui(2)=particlei%uy
  ui(3)=particlei%uz
  uj(1)=particlej%ux
  uj(2)=particlej%uy
  uj(3)=particlej%uz
  if (particlei%rod .and. particlej%rod) then
    call potential(ui, uj, rij, potE, overlap)
  else if (particlei%rod) then
    call gblj_potential(ui, rij, potE, overlap)
  else if (particlej%rod) then
    call gblj_potential(uj, -rij, potE, overlap)
  else
    potE = lj_potential(sqrt(dot_product(rij, rij)))
    overlap = .false.
  end if
end subroutine


pure function pair_force(particlei, particlej, rij)
  type(particledat), intent(in) :: particlei, particlej
  real(dp), intent(in) :: rij(3)
  real(dp) :: pair_force(3)
  real(dp), dimension(3) :: ui, uj
  ui(1)=particlei%ux
  ui(2)=particlei%uy
  ui(3)=particlei%uz
  uj(1)=particlej%ux
  uj(2)=particlej%uy
  uj(3)=particlej%uz
  if (particlei%rod .and. particlej%rod) then
    pair_force = gb_force(ui, uj, rij)
  else if (particlei%rod) then
    pair_force = gblj_force(ui, rij)
  else if (particlej%rod) then
    pair_force = gblj_force(uj, -rij)
  else
    pair_force = lj_force(rij)
  end if
  
end function 

end module
