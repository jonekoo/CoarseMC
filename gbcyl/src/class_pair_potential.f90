!> Module responsible for calculations of pair interactions between
!! particles.
module class_pair_potential
  use iso_fortran_env
  use num_kind
  use m_sphere_interaction
  use m_rod_interaction
  use m_rodsphere_potential
  use m_gayberne
  use m_gblj
  use m_lj
  use particle, only: particledat, pair_interaction, position
  use class_poly_box, only: poly_box, minimage
  use class_parameterizer
  use class_parameter_writer
  implicit none
  
  type, extends(pair_interaction) :: conditional_pair_interaction
     real(dp) :: cutoff = 5.5_dp
     class(rod_interaction), allocatable :: p_rod_ia
     class(sphere_interaction), allocatable :: p_sphere_ia
     class(rodsphere_potential), allocatable :: p_rodsphere
   contains
     procedure :: pair_potential => cpi_pair_potential
     procedure :: pair_force => cpi_pair_force
     procedure :: get_cutoff => cpi_get_cutoff
  end type conditional_pair_interaction

  type(conditional_pair_interaction), save, private :: prototype
  logical, save, private :: is_initialized = .false.

  interface assignment(=)
     module procedure assign
  end interface assignment(=)
  
contains

  function create_conditional_interaction() result(res)
    type(conditional_pair_interaction) :: res
    if (is_initialized) then
       res = prototype
    else
       stop 'class_pair_potential not initialized!'
    end if
  end function create_conditional_interaction
  
  subroutine assign(lhs, rhs)
    type(conditional_pair_interaction), intent(out) :: lhs
    type(conditional_pair_interaction), intent(in) :: rhs
    allocate(lhs%p_rod_ia, source=rhs%p_rod_ia)
    allocate(lhs%p_sphere_ia, source=rhs%p_sphere_ia)
    allocate(lhs%p_rodsphere, source=rhs%p_rodsphere)
  end subroutine assign

  pure function cpi_get_cutoff(this) result(cutoff)
    class(conditional_pair_interaction), intent(in) :: this
    real(dp) :: cutoff
    cutoff = this%cutoff
  end function cpi_get_cutoff
  
  !> Initializes the dependencies of this module.
  !! 
  !! @param[in] reader the object responsible for reading the parameters.
  !! 
  subroutine pp_init(reader)
    type(parameterizer), intent(in) :: reader
    !allocate(gayberne :: p_rod_ia)
    character(len=20) :: rod_type_str, sphere_type_str, rodsphere_type_str
    rod_type_str = 'gayberne'
    call getparameter(reader, 'r_cutoff', prototype%cutoff)
    call getparameter(reader, 'rod_potential', rod_type_str)
    if (trim(adjustl(rod_type_str)) == 'gayberne') then
       allocate(prototype%p_rod_ia, source=gayberne(reader))
    else
       write(error_unit, *) 'Error: unknown rod_potential', rod_type_str
       stop 1
    end if

    sphere_type_str = 'lj'
    call getparameter(reader, 'sphere_potential', sphere_type_str)
    if (trim(adjustl(sphere_type_str)) == 'lj') then
       allocate(prototype%p_sphere_ia, source=lj_potential(reader))
    else
       write(error_unit, *) 'Error: unknown sphere_potential', sphere_type_str
       stop 1
    end if

    rodsphere_type_str = 'gblj'
    call getparameter(reader, 'rodsphere_potential', rodsphere_type_str)
    if(trim(adjustl(rodsphere_type_str)) == 'gblj') then
       allocate(prototype%p_rodsphere, source=gblj_potential(reader))
    else
       write(error_unit, *) 'Error: unknown rodsphere_potential', &
            rodsphere_type_str
       stop 1
    end if

    is_initialized = .true.
  end subroutine pp_init
  
  !> Hands the parameter @p writer over to the dependencies of this
  !! module.
  !! 
  !! @param[in] writer the object responsible for formatting the output and
  !! handling the file. 
  !!
  subroutine pp_writeparameters(writer)
    type(parameter_writer), intent(inout) :: writer
    call writeparameter(writer, 'r_cutoff', prototype%cutoff)
    select type (p => prototype%p_rod_ia)
    type is (gayberne)
       call writeparameter(writer, 'rod_potential', 'gayberne')
    class default
       call writecomment(writer, 'Warning: rod_potential not recognized.')
    end select
    call prototype%p_rod_ia%writeparameters(writer)

    select type (p => prototype%p_sphere_ia)
    type is (lj_potential)
       call writeparameter(writer, 'sphere_potential', 'lj')
    class default
       call writecomment(writer, 'Warning: sphere_potential not recognized.')
    end select
    call prototype%p_sphere_ia%writeparameters(writer)

    select type (p => prototype%p_rodsphere)
    type is (gblj_potential)
       call writeparameter(writer, 'rodsphere_potential', 'gblj')
    class default
       call writecomment(writer, 'Warning: rodsphere_potential not recognized.')
    end select
    call prototype%p_rodsphere%writeparameters(writer)
  end subroutine pp_writeparameters
  

  pure subroutine cpi_pair_potential_2(this, particlei, particlej, &
       simbox, energy, err)
    class(conditional_pair_interaction), intent(in) :: this
    type(particledat), intent(in) :: particlei, particlej
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    err = 0
    call this%pair_potential(particlei, particlej, &
         minimage(simbox, position(particlej) - position(particlei)), energy, &
         err)
  end subroutine cpi_pair_potential_2
  
  !> Calculates the interaction energy of a pair of particles. 
  !! 
  !! @param[in] particlei,particlej the pair of particles.
  !! @param[in] rij is the (minimum image) vector from particle i to particle j.
  !! @param[out] potE the interaction energy.
  !! @param[out] overlap is true if the two particles overlap each other.
  !! 
  pure subroutine cpi_pair_potential(this, particlei, particlej, rij, energy,&
       err)
    class(conditional_pair_interaction), intent(in) :: this
    type(particledat), intent(in) :: particlei 
    type(particledat), intent(in) :: particlej
    real(dp), intent(in) :: rij(3)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    real(dp), dimension(3) :: ui, uj
    ui(1) = particlei%ux
    ui(2) = particlei%uy
    ui(3) = particlei%uz
    uj(1) = particlej%ux
    uj(2) = particlej%uy
    uj(3) = particlej%uz
    err = 0
    if (particlei%rod .and. particlej%rod) then
       call this%p_rod_ia%potential(ui, uj, rij, energy, overlap)
    else if (particlei%rod) then
       call this%p_rodsphere%potential(ui, rij, energy, overlap)
    else if (particlej%rod) then
       call this%p_rodsphere%potential(uj, -rij, energy, overlap)
    else
       call this%p_sphere_ia%potential(sqrt(dot_product(rij, rij)), &
            energy, overlap)
    end if
    if (overlap) err = 1
  end subroutine cpi_pair_potential
  

  !> The force acting on @p particlej caused by @p particlei.
  !!
  !! @param[in] particlei, particlej the pair of particles.
  !! @param[in] rij = rj - ri, vector from particle i to particle j.
  !!
  !! @return the force acting on @p particlej
  !!
  pure function cpi_pair_force(this, particlei, particlej, rij) result(f)
    class(conditional_pair_interaction), intent(in) :: this
    type(particledat), intent(in) :: particlei, particlej
    real(dp), intent(in) :: rij(3)
    real(dp) :: f(3)
    real(dp), dimension(3) :: ui, uj
    ui(1) = particlei%ux
    ui(2) = particlei%uy
    ui(3) = particlei%uz
    uj(1) = particlej%ux
    uj(2) = particlej%uy
    uj(3) = particlej%uz
    if (particlei%rod .and. particlej%rod) then
       f = this%p_rod_ia%force(ui, uj, rij)
    else if (particlei%rod) then
       f = this%p_rodsphere%force(ui, rij)
    else if (particlej%rod) then
       f = this%p_rodsphere%force(uj, -rij)
    else
       f = this%p_sphere_ia%force(rij)
    end if
  end function cpi_pair_force

end module class_pair_potential
