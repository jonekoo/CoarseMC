!> Module responsible for calculations of pair interactions between
!! particles.
module class_pair_potential
  use iso_fortran_env
  use nrtype
  use m_sphere_interaction
  use m_rod_interaction
  use m_rodsphere_potential
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

  class(rod_interaction), allocatable :: p_rod_ia
  class(sphere_interaction), allocatable :: p_sphere_ia
  class(rodsphere_potential), allocatable :: p_rodsphere
  
contains
  
  !> Initializes the dependencies of this module.
  !! 
  !! @param[in] reader the object responsible for reading the parameters.
  !! 
  subroutine pp_init(reader)
    type(parameterizer), intent(in) :: reader
    !allocate(gayberne :: p_rod_ia)
    character(len=20) :: rod_type_str, sphere_type_str, rodsphere_type_str
    rod_type_str = 'gayberne'
    call getparameter(reader, 'rod_potential', rod_type_str)
    if (trim(adjustl(rod_type_str)) == 'gayberne') then
       allocate(p_rod_ia, source=gayberne(reader))
    else
       write(error_unit, *) 'Error: unknown rod_potential', rod_type_str
       stop 1
    end if

    sphere_type_str = 'lj'
    call getparameter(reader, 'sphere_potential', sphere_type_str)
    if (trim(adjustl(sphere_type_str)) == 'lj') then
       allocate(p_sphere_ia, source=lj_potential(reader))
    else
       write(error_unit, *) 'Error: unknown sphere_potential', sphere_type_str
       stop 1
    end if

    rodsphere_type_str = 'gblj'
    call getparameter(reader, 'rodsphere_potential', rodsphere_type_str)
    if(trim(adjustl(rodsphere_type_str)) == 'gblj') then
       allocate(p_rodsphere, source=gblj_potential(reader))
    else
       write(error_unit, *) 'Error: unknown rodsphere_potential', &
            rodsphere_type_str
       stop 1
    end if
    
  end subroutine pp_init
  
  !> Hands the parameter @p writer over to the dependencies of this
  !! module.
  !! 
  !! @param[in] writer the object responsible for formatting the output and
  !! handling the file. 
  !!
  subroutine pp_writeparameters(writer)
    type(parameter_writer), intent(inout) :: writer
    select type (p_rod_ia)
    type is (gayberne)
       call writeparameter(writer, 'rod_potential', 'gayberne')
    class default
       call writecomment(writer, 'Warning: rod_potential not recognized.')
    end select
    call p_rod_ia%writeparameters(writer)

    select type (p_sphere_ia)
    type is (lj_potential)
       call writeparameter(writer, 'sphere_potential', 'lj')
    class default
       call writecomment(writer, 'Warning: sphere_potential not recognized.')
    end select
    call p_sphere_ia%writeparameters(writer)

    select type (p_rodsphere)
    type is (gblj_potential)
       call writeparameter(writer, 'rodsphere_potential', 'gblj')
    class default
       call writecomment(writer, 'Warning: rodsphere_potential not recognized.')
    end select
    call p_rodsphere%writeparameters(writer)
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
       call p_rod_ia%potential(ui, uj, rij, potE, overlap)
    else if (particlei%rod) then
       call p_rodsphere%potential(ui, rij, potE, overlap)
    else if (particlej%rod) then
       call p_rodsphere%potential(uj, -rij, potE, overlap)
    else
       call p_sphere_ia%potential(sqrt(dot_product(rij, rij)), potE, overlap)
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
       pair_force = p_rod_ia%force(ui, uj, rij)
    else if (particlei%rod) then
       pair_force = p_rodsphere%force(ui, rij)
    else if (particlej%rod) then
       pair_force = p_rodsphere%force(uj, -rij)
    else
       pair_force = p_sphere_ia%force(rij)
    end if
  end function pair_force

end module class_pair_potential
