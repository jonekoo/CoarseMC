!> Module responsible for calculations of pair interactions between
!! particles.
module class_pair_potential
  use nrtype
  use m_gayberne
  use m_gblj
  use m_lj
  use particle
  use class_parameterizer
  use class_parameter_writer
  implicit none
  private
  
  public :: pair_potential
  public :: pp_init
  public :: pp_writeparameters
  public :: pair_force
  
  type(gayberne) :: gb
  type(gblj_potential) :: gblj
  type(lj_potential) :: lj
  
contains
  
  !> Initializes the dependencies of this module.
  !! 
  !! @param[in] reader the object responsible for reading the parameters.
  !! 
  subroutine pp_init(reader)
    type(parameterizer), intent(in) :: reader
    gb = gayberne(reader)
    lj = lj_potential(reader)
    gblj = gblj_potential(reader)
  end subroutine pp_init
  
  !> Hands the parameter @p writer over to the dependencies of this
  !! module.
  !! 
  !! @param[in] writer the object responsible for formatting the output and
  !! handling the file. 
  !!
  subroutine pp_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call gb%writeparameters(writer)
    call lj%writeparameters(writer)
    call gblj%writeparameters(writer)
  end subroutine pp_writeparameters
  
  
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
       call gb%potential(ui, uj, rij, potE, overlap)
    else if (particlei%rod) then
       call gblj%potential(ui, rij, potE, overlap)
    else if (particlej%rod) then
       call gblj%potential(uj, -rij, potE, overlap)
    else
       potE = lj%potential(sqrt(dot_product(rij, rij)))
       overlap = .false.
    end if
  end subroutine pair_potential
  

  !> The force acting on @p particlej caused by @p particlei.
  !!
  !! @param[in] particlei, particlej the pair of particles.
  !! @param[in] rij = rj - ri, vector from particle i to particle j.
  !!
  !! @return the force acting on @p particlej
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
       pair_force = gb%force(ui, uj, rij)
    else if (particlei%rod) then
       pair_force = gblj%force(ui, rij)
    else if (particlej%rod) then
       pair_force = gblj%force(uj, -rij)
    else
       pair_force = lj%force(rij)
    end if
  end function pair_force

end module class_pair_potential
