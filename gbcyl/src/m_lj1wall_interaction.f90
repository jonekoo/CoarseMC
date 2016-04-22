module m_lj1wall_interaction
  use m_particle, only: point, single_interaction 
  use num_kind
  use class_poly_box
  use ljcylinder
  use json_module
  use m_json_wrapper, only: get_parameter
  implicit none

  !> Defines the interaction between a Lennard-Jones particle and the
  !! wall of a cylindrical cavity confining the particle. The wall
  !! consists of smoothly and evenly distributed "virtual" LJ particles.
  !!
  !! @see X. Zhang et al., Fluid Phase Equilibria 218(2), 239-246, 2004.
  !!
  type, extends(single_interaction) :: lj1wall_interaction
     !> The relative strength of the attractive term as compared to the
     !! repulsive term: U = U_repulsive + alpha * U_attractive
     real(dp) :: alpha = 1.

     !> The well-depth parameter of the interaction between the wall and the
     !! LJ site.
     real(dp) :: eps = 1.
     
     !> The range parameter of the interaction between the LJ particle and the
     !! wall.
     real(dp) :: sig = 1.
     
     !> Number density of virtual LJ particles in the wall.
     real(dp) :: wall_density = 1.
   contains

     !> Computes the potential energy of the interaction between
     !! a LJ particle and the wall.
     procedure :: potential => lj1wall

     !> Computes the force acting on the given LJ particle by the wall.
     procedure :: force => lj1wall_force
     
     !> Serializes the interaction to JSON.
     procedure :: to_json => lj1wall_ia_to_json

     !> Produces a sample of the potential energies with this
     !! interaction.
     procedure :: sample => lj1wall_sample
  end type lj1wall_interaction


  interface lj1wall_interaction
     module procedure lj1wall_interaction_from_json
  end interface lj1wall_interaction

  
contains

  !> Creates the lj1wall_interaction using parameters read from JSON
  !! contained in the @p json_val.
  function lj1wall_interaction_from_json(json_val) result(res)
    type(json_value), pointer, intent(in) :: json_val
    type(lj1wall_interaction) :: res
    call get_parameter(json_val, 'alpha', res%alpha, error_lb=0._dp) 
    call get_parameter(json_val, 'sig', res%sig, error_lb=0._dp)
    call get_parameter(json_val, 'wall_density', res%wall_density, &
         error_lb=0._dp)
    call get_parameter(json_val, 'eps', res%eps, error_lb=0._dp)
  end function


  !> Writes the parameters of @p this interaction to the JSON value
  !! @p json_val.
  subroutine lj1wall_ia_to_json(this, json_val)
    class(lj1wall_interaction), intent(in) :: this
    type(json_value), intent(inout), pointer :: json_val    
    call json_add(json_val, 'type', 'lj1wall_interaction')
    call json_add(json_val, 'wall_density', this%wall_density)
    call json_add(json_val, 'sig', this%sig)
    call json_add(json_val, 'eps', this%eps)
    call json_add(json_val, 'alpha', this%alpha)
  end subroutine


  !> Calculates the interaction energy of a Lennard-Jones (LJ) particle
  !! and the wall of a cylindrical cavity.
  !!
  !! @param this the LJ-wall interaction.
  !! @param particlei the LJ particle.
  !! @param simbox the simulation box defining the radius of the
  !!        cavity.
  !! @param energy the interaction energy.
  !! @param err == 1 if the @p particlei has penetrated the
  !! wall too much.
  !! 
  pure subroutine lj1wall(this, particlei, simbox, energy, err)
    class(lj1wall_interaction), intent(in) :: this
    class(point), intent(in) :: particlei
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    real(dp) :: r
    real(dp) :: r_cylinder
    err = 0
    energy = 0._dp
    r_cylinder = getx(simbox)/2._dp 
    r = sqrt(particlei%x**2 + particlei%y**2)
    if(r >= r_cylinder) then
      err = 1
      return
    end if
    energy = ljcylinderpotential(this%eps, this%wall_density, &
         this%sig, this%alpha, r, r_cylinder)
  end subroutine


  !> Returns the force exerted on a LJ particle @p particlei by the
  !! wall of the cylindrical cavity due to @p this interaction. The
  !! radius of the cavity is defined by the @p simbox.
  pure function lj1wall_force(this, particlei, simbox) result(f)
    class(lj1wall_interaction), intent(in) :: this
    class(point), intent(in) :: particlei
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)

    real(dp) :: r
    real(dp) :: r_cylinder
    r_cylinder = getx(simbox)/2._dp 
    r = sqrt(particlei%x**2 + particlei%y**2)
    !! Set the direction of f:
    f = [particlei%x, particlei%y, 0._dp] / r
    !! What should be done when r is near zero in the division above?
    f = f * ljcylinderforce(this%eps, this%wall_density, &
         this%sig, this%alpha, r, r_cylinder)
  end function


  !> Creates a sample of the possible energies of @p this interaction
  !! with wall-particle distances @p r and stores them to the JSON
  !! value @p json_val.
  !!
  subroutine lj1wall_sample(this, json_val, r, simbox)
    class(lj1wall_interaction), intent(in) :: this
    type(json_value), pointer :: json_val
    real(dp), intent(in) :: r(:)
    type(poly_box), intent(in) :: simbox
    type(point) :: lj
    type(json_value), pointer :: child
    real(dp) :: energy
    integer :: err, i, k
    type(json_value), pointer :: r_json, energy_json

    ! LJ-wall interaction
    call json_create_object(child, 'LJ' // '-LJwall')
    call json_create_array(r_json, 'x')
    call json_create_array(energy_json, 'energy')
    do k = 1, size(r)
       call lj%set_position([simbox%lx / 2 - r(k), 0._dp, 0._dp])
       call this%potential(lj, simbox, energy, err)
       if (err == 0) then
          call json_add(r_json, '', r(k))
          call json_add(energy_json, '', energy)
       end if
    end do
    call json_add(child, r_json)
    call json_add(child, energy_json)
    call json_add(json_val, child)
  end subroutine lj1wall_sample
  
  
end module
