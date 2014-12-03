!> Module responsible for calculations of pair interactions between
!! particles.
module class_pair_potential
use nrtype
use gayberne
use lj, only: lj_potential, lj_init, lj_writeparameters, lj_force
use gblj, only: gblj_potential, gblj_init, gblj_force, gblj_writeparameters
use particle
use class_parameterizer
use class_parameter_writer
implicit none
private

public :: pair_potential
public :: pp_init
public :: pp_writeparameters
public :: pair_force

contains

!> Initializes the dependencies of this module.
!! 
!! @param[in] reader the object responsible for reading the parameters.
!! 
subroutine pp_init(reader)
  type(parameterizer), intent(in) :: reader
  call gayberne_init(reader)
  call lj_init(reader)
  call gblj_init(reader)
end subroutine

!> Hands the parameter @p writer over to the dependencies of this
!! module.
!! 
!! @param[in] writer the object responsible for formatting the output and
!! handling the file. 
!!
subroutine pp_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call gb_writeparameters(writer)
  call lj_writeparameters(writer)
  call gblj_writeparameters(writer)
end subroutine


!> Calculates the interaction energy of a pair of particles. 
!! 
!! @param[in] particlei,particlej the pair of particles.
!! @param[in] rij is the (minimum image) vector from particle i to particle j.
!! @param[out] potE the interaction energy.
!! @param[out] overlap is true if the two particles overlap each other.
!! 
pure subroutine pair_potential(particlei, particlej, rij, potE, overlap)
  type(particledat), intent(in) :: particlei 
  type(particledat), intent(in) :: particlej
  real(dp), intent(in) :: rij(3)
  real(dp), intent(out) :: potE
  logical, intent(out) :: overlap
  real(dp), dimension(3) :: ui, uj
  ui(1) = particlei%ux
  ui(2) = particlei%uy
  ui(3) = particlei%uz
  uj(1) = particlej%ux
  uj(2) = particlej%uy
  uj(3) = particlej%uz
  if (particlei%rod .and. particlej%rod) then
    call gb_potential(ui, uj, rij, potE, overlap)
  else if (particlei%rod) then
    call gblj_potential(ui, rij, potE, overlap)
  else if (particlej%rod) then
    call gblj_potential(uj, -rij, potE, overlap)
  else
    potE = lj_potential(sqrt(dot_product(rij, rij)))
    overlap = .false.
  end if
end subroutine


!> The force acting on @p particlei caused by @p particlej.
!!
!! @param[in] particlei,particlej the pair of particles.
!! @param[in] rij the minimum image distance of the two particles.
!!
!! @return the force acting on @p particlei
!!
pure function pair_force(particlei, particlej, rij)
  type(particledat), intent(in) :: particlei, particlej
  real(dp), intent(in) :: rij(3)
  real(dp) :: pair_force(3)
  real(dp), dimension(3) :: ui, uj
  ui(1) = particlei%ux
  ui(2) = particlei%uy
  ui(3) = particlei%uz
  uj(1) = particlej%ux
  uj(2) = particlej%uy
  uj(3) = particlej%uz
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
