!> Module for computing particle-wall interactions for a Lennard-Jones
!! (LJ)  or Gay-Berne (GB) particle in a cylindrical cavity. The cavity
!! walls are smooth and consist of evenly distributed LJ particles. A
!! GB particle interacts with the wall through two embedded LJ sites in
!! the particle.
!! 
!! @see Micheletti et. al. Journal of Chemical Physics 123, 224705.
!!
module particlewall
  use m_particle, only: particle, single_interaction 
  use num_kind
  use class_poly_box
  use class_parameterizer
  use class_parameter_writer
  use ljcylinder
  use json_module
  use m_json_wrapper, only: get_parameter
  use m_rod, only: rod
  use m_point, only: point
  implicit none
    
  type, extends(single_interaction) :: ljwall_interaction
     !> The relative strength of attractive term as compared to the
     !! repulsive term for the LJ site A of the GB particle. 
     real(dp) :: alpha_a = 1.
     
     !> The relative strength of attractive term as compared to the
     !! repulsive term for the LJ site B of the GB particle. 
     real(dp) :: alpha_b = 1.
     
     !> The well-depth parameter of the interaction between the wall and the
     !! LJ site.
     real(dp) :: eps = 1.
     
     !> The range parameter of the interaction between the LJ site and the
     !! wall.
     real(dp) :: sig = 1.
     
     !> The distance of LJ interaction sites from the Gay-Berne particle
     !! center along the unique axis.
     real(dp) :: LJdist = 1.7 
  
     !> If true, the wall interaction favors uniform alignment of GB 
     !! particles along the z-axis of the simulation box.
     logical :: isuniformalignment = .false.
     
     !> Well-depth parameter for the interaction of a LJ particle with the
     !! wall.
     real(dp) :: epswall_lj = 1.
     
     !> The strength of the attractive term with respect to the repulsive
     !! term for the interaction of a LJ particle with the wall.
     real(dp) :: alpha_lj = 1.
     
     !> Range parameter for the interaction of a LJ particle with the
     !! wall.
     real(dp) :: sigwall_lj = 0.8
     
     !> Number density of virtual LJ particles in the wall.
     real(dp) :: wall_density = 1.
   contains

     !> Computes the potential energy of the interaction between
     !! a particle and the wall.
     procedure :: potential => ljwall_ia_potential

     !> Computes the force acting on the given particle by the wall.
     procedure :: force => ljwall_ia_force

     !> Serializes the interaction to JSON.
     procedure :: to_json => ljwall_ia_to_json

     !> Computes the interaction energy of a rod with the wall.
     procedure :: rod_potential => gbwall

     !> Computes the interaction energy of a point with the wall.
     procedure :: point_potential => ljwall

     !> Returns the force acting on a rod by the wall.
     procedure :: rod_force => gbwall_force

     !> Returns the force acting on a point by the wall.
     procedure :: point_force => ljwall_force

     !> Returns the distances of the interactions sites in a rod
     !! from the cylinder axis.
     procedure :: rarb

     !> Produces a sample of the potential energies with this
     !! interaction.
     procedure :: sample => ljwall_sample
 end type ljwall_interaction

 !> Constructors for the ljwall_interaction
 interface ljwall_interaction
    module procedure ljwall_interaction_from_json, ljwall_interaction_w_reader
 end interface ljwall_interaction
 
contains 

  !> Creates a ljwall_interaction by reading its parameters with
  !! @p reader.
  function ljwall_interaction_w_reader(reader) result(res)
    type(ljwall_interaction) :: res
    type(parameterizer), intent(in) :: reader
    real(dp) :: Kw_LJ = 5.48819_dp, Kw = 8._dp
    logical :: found
    call getparameter(reader, 'alpha_A', res%alpha_a) 
    call getparameter(reader, 'alpha_B', res%alpha_b) 
    call getparameter(reader, 'LJ_dist', res%LJdist) 
    call getparameter(reader, 'is_uniform_alignment', res%isuniformalignment)
    call getparameter(reader, 'sigwall', res%sig)
    call getparameter(reader, 'alpha_LJ', res%alpha_lj)
    call getparameter(reader, 'sigwall_LJ', res%sigwall_lj)
    call getparameter(reader, 'wall_density', res%wall_density, found)
    if (.not. found) then
       !! Try to read old parameters
       call getparameter(reader, 'Kw', Kw)
       res%eps = Kw / 8._dp
       call getparameter(reader, 'Kw_LJ', Kw_LJ)
       res%epswall_lj = Kw_LJ / (Kw * (res%sigwall_lj/res%sig)**3) 
    else
       !! Use the new parameters epswall_lj and epswall
       call getparameter(reader, 'epswall_LJ', res%epswall_lj)
       call getparameter(reader, 'epswall', res%eps)
    end if
  end function



  !> Creates the ljwall_interaction using parameters read from JSON
  !! contained in the @p json_val.
  function ljwall_interaction_from_json(json_val) result(res)
    type(json_value), pointer, intent(in) :: json_val
    type(ljwall_interaction) :: res
    call get_parameter(json_val, 'alpha_A', res%alpha_a, error_lb=0._dp) 
    call get_parameter(json_val, 'alpha_B', res%alpha_b, error_lb=0._dp) 
    call get_parameter(json_val, 'LJ_dist', res%LJdist, error_lb=0._dp) 
    call get_parameter(json_val, 'is_uniform_alignment', res%isuniformalignment)
    call get_parameter(json_val, 'sigwall', res%sig, error_lb=0._dp)
    call get_parameter(json_val, 'alpha_LJ', res%alpha_lj, error_lb=0._dp)
    call get_parameter(json_val, 'sigwall_LJ', res%sigwall_lj, error_lb=0._dp)
    call get_parameter(json_val, 'wall_density', res%wall_density, &
         error_lb=0._dp)
    call get_parameter(json_val, 'epswall_LJ', res%epswall_lj, error_lb=0._dp)
    call get_parameter(json_val, 'epswall', res%eps, error_lb=0._dp)
  end function


  !> Writes the parameters of this module with the format and output
  !! unit defined by @p writer.
  subroutine ljwall_ia_to_json(this, json_val)
    class(ljwall_interaction), intent(in) :: this
    type(json_value), intent(inout), pointer :: json_val
    
    call json_add(json_val, 'type', 'ljwall_interaction')
    call json_add(json_val, 'wall_density', this%wall_density)

    call json_add(json_val, 'rodwall', 'ljdimer-wall')
    call json_add(json_val, 'alpha_A', this%alpha_a) 
    call json_add(json_val, 'alpha_B', this%alpha_b) 
    call json_add(json_val, 'LJ_dist', this%LJdist) 
    call json_add(json_val, 'is_uniform_alignment', this%isuniformalignment)
    call json_add(json_val, 'sigwall', this%sig)
    call json_add(json_val, 'epswall', this%eps)

    call json_add(json_val, 'pointwall', 'ljcylinder')
    call json_add(json_val, 'alpha_LJ', this%alpha_lj)
    call json_add(json_val, 'sigwall_LJ', this%sigwall_lj)
    call json_add(json_val, 'epswall_LJ', this%epswall_lj)
  end subroutine


  !> Calculates the potential @p energy for a @p particle inside a
  !! cylindrical Lennard-Jones cavity. The radius of the cavity is
  !! defined by @p simbox. If the particle is too close to the wall or
  !! inside the wall @p overlap == true.
  pure subroutine ljwall_ia_potential(this, particlei, simbox, energy, err)
    class(ljwall_interaction), intent(in) :: this
    class(particle), intent(in) :: particlei
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    err = 0
    select type (particlei)
    type is (rod)
      call this%rod_potential(particlei, simbox, energy, overlap)
    type is (point)
      call this%point_potential(particlei, simbox, energy, overlap)
    class default
      err = 66
   end select
   if (overlap) err = 1
 end subroutine


  !> Returns the force acting on @p particlei by the wall of a
  !! cylindrical cavity. Radius of the cavity is defined by @p simbox.
  !! See module description for details.
  pure function ljwall_ia_force(this, particlei, simbox) result(f)
    class(ljwall_interaction), intent(in) :: this
    class(particle), intent(in) :: particlei
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)
    f = 0.
    select type (particlei)
    type is (rod)
      f = this%rod_force(particlei, simbox)
    type is (point)
      f = this%point_force(particlei, simbox)
    end select
  end function


  !> Calculates the interaction energy of a Lennard-Jones (LJ) particle
  !! and the wall of a cylindrical cavity.
  !! 
  !! @param ljparticle the LJ particle.
  !! @param simbox the simulation box defining the radius of the
  !!        cavity.
  !! @param energy the interaction energy.
  !! @param overlap is true if the @p ljparticle has penetrated the
  !! wall too much.
  !! 
  pure subroutine ljwall(this, ljparticle, simbox, energy, overlap)
    class(ljwall_interaction), intent(in) :: this
    class(point), intent(in) :: ljparticle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    real(dp) :: r
    real(dp) :: r_cylinder
    overlap = .false.
    energy = 0._dp
    r_cylinder = getx(simbox)/2._dp 
    r = sqrt(ljparticle%x**2 + ljparticle%y**2)
    if(r >= r_cylinder) then
      overlap = .true.
      return
    end if
    energy = ljcylinderpotential(this%epswall_lj, this%wall_density, &
         this%sigwall_lj, this%alpha_lj, r, r_cylinder)
  end subroutine


  !> Returns the force exerted on a LJ particle @p ljparticle by the
  !! wall of the cylindrical cavity. The radius of the cavity is
  !! defined by the @p simbox.
  pure function ljwall_force(this, ljparticle, simbox) result(f)
    class(ljwall_interaction), intent(in) :: this
    class(point), intent(in) :: ljparticle
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)

    real(dp) :: r
    real(dp) :: r_cylinder
    r_cylinder = getx(simbox)/2._dp 
    r = sqrt(ljparticle%x**2 + ljparticle%y**2)
    !! Set the direction of f:
    f = [ljparticle%x, ljparticle%y, 0._dp] / r
    !! What should be done when r is near zero in the division above?
    f = f * ljcylinderforce(this%epswall_lj, this%wall_density, &
         this%sigwall_lj, this%alpha_lj, r, r_cylinder)
  end function


  !> Calculates the potential energy of a rodlike particle which
  !! interacts with the wall of a cylindrical cavity. The particle
  !! interacts with the wall via two embedded Lennard-Jones interaction
  !! sites.
  !! 
  !! @param gbparticle the rodlike particle.
  !! @param simbox the simulation box defining the dimensions of the
  !!        cavity.
  !! @param energy the interaction energy of the wall and the particle.
  !! @param overlap is true if the @p gbparticle has penetrated the
  !!        wall too much.
  !! 
  !! @see D. Micheletti et al. J. Chem. Phys. 123, 224705, 2005 for the
  !! interaction site model.
  !!
  pure subroutine gbwall(this, gbparticle, simbox, energy, overlap)
    implicit none
    class(ljwall_interaction), intent(in) :: this
    class(rod), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    real(dp) :: rsite_a, rsite_b
    real(dp) :: fu, r_cylinder
    overlap = .false.
    energy = 0._dp
    r_cylinder = getx(simbox)/2._dp 
    call this%rarb(gbparticle, rsite_a, rsite_b)
    if(rsite_a >= r_cylinder .or. rsite_b >= r_cylinder) then
      overlap = .true.
      return
    end if
    if (this%isuniformalignment) then
      fu = angular(gbparticle)  
    else 
      fu = 1._dp
    end if
    energy = fu * (ljcylinderpotential(this%eps, this%wall_density, this%sig,&
         this%alpha_a, rsite_a, r_cylinder) + &
         ljcylinderpotential(this%eps, this%wall_density, this%sig, &
         this%alpha_b, rsite_b, r_cylinder))
  end subroutine gbwall


  !> Returns the force exerted on a GB particle @p gbparticle by the
  !! wall of a cylindrical cavity. The @p gbparticle interacts with the
  !! wall via two embedded LJ sites.
  pure function gbwall_force(this, gbparticle, simbox) result(f)
    class(ljwall_interaction), intent(in) :: this
    class(rod), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)
    real(dp) :: f_a(3), f_b(3)
    real(dp) :: rsite_a, rsite_b
    real(dp) :: fu, r_cylinder
    r_cylinder = getx(simbox) / 2. 
    call this%rarb(gbparticle, rsite_a, rsite_b)
    if (this%isuniformalignment) then
      fu = angular(gbparticle)  
    else 
      fu = 1.
    end if

    f_a = gbparticle%position() + gbparticle%orientation() * this%LJdist
    f_a(3) = 0._dp
    f_a = f_a / sqrt(dot_product(f_a, f_a))

    f_b = gbparticle%position() - gbparticle%orientation() * this%LJdist
    f_b(3) = 0._dp
    f_b = f_b / sqrt(dot_product(f_b, f_b))

    f = fu * (f_a * ljcylinderforce(this%eps, this%wall_density, this%sig,&
         this%alpha_a, rsite_a, r_cylinder) + &
         f_b * ljcylinderforce(this%eps, this%wall_density, this%sig,&
         this%alpha_b, rsite_b, r_cylinder))
  end function


  !> Returns the distances from the cavity axis @p ra and @p ra for the
  !! interaction sites A and B, respectively. The interaction sites are
  !! embedded in @p particle.
  pure subroutine rarb(this, particle, ra, rb)
    class(ljwall_interaction), intent(in) :: this
    class(rod), intent(in) :: particle
    real(dp), intent(out) :: ra, rb
    real(dp) :: xa, ya, xb, yb
    xa = particle%x + this%LJdist * particle%ux
    ya = particle%y + this%LJdist * particle%uy
    xb = particle%x - this%LJdist * particle%ux
    yb = particle%y - this%LJdist * particle%uy
    ra = sqrt(xa**2 + ya**2)
    rb = sqrt(xb**2 + yb**2)
  end subroutine


  !> Returns the angular dependence of the potential with respect to the 
  !! cylinder axis (z-direction). 
  !!
  !! @param particle the particle to which the potential is calculated. 
  !! 
  pure real(dp) function angular(particle)
    class(rod), intent(in) :: particle
    angular = (particle%uz)**2
  end function angular


  !> Creates a sample of the possible energies at distances @p r with
  !! rod and point particles and stores them to the JSON value
  !! @p json_val. For the rod, three different orientations along x, y
  !! and z axes are used.
  !!
  subroutine ljwall_sample(this, json_val, r, simbox)
    class(ljwall_interaction), intent(in) :: this
    type(json_value), pointer :: json_val
    real(dp), intent(in) :: r(:)
    type(poly_box), intent(in) :: simbox
    type(rod) :: rodsamples(3)
    type(point) :: lj
    type(json_value), pointer :: child
    real(dp) :: energy
    integer :: err, i, j, k
    type(json_value), pointer :: r_json, energy_json
    character(len=3, kind=CK), parameter :: descriptions(3) = &
         ['GBx', 'GBy', 'GBz']

    ! First the GB-Wall interactions
    call rodsamples(1)%set_orientation([1._dp, 0._dp, 0._dp])
    call rodsamples(2)%set_orientation([0._dp, 1._dp, 0._dp])
    call rodsamples(3)%set_orientation([0._dp, 0._dp, 1._dp])
    call json_create_object(json_val, 'ljwall_interaction')
    do i = 1, size(rodsamples)   
       call json_create_object(child, descriptions(i) // '-LJwall')
       call json_create_array(r_json, 'x')
       call json_create_array(energy_json, 'energy')
       do k = 1, size(r)
          call rodsamples(i)%set_position([simbox%lx / 2 - r(k), 0._dp, 0._dp])
          call this%potential(rodsamples(i), &
               simbox, energy, err)
          if (err == 0) then
             call json_add(r_json, '', r(k))
             call json_add(energy_json, '', energy)
          end if
       end do
       call json_add(child, r_json)
       call json_add(child, energy_json)
       call json_add(json_val, child)
    end do
    
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
  end subroutine ljwall_sample
  
end module particlewall
