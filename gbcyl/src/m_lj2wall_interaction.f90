!> Particle-wall interactions for Lennard-Jones and rod-like particles
!! in a cylindrical cavity. The cavity walls are smooth and consist of
!! evenly distributed LJ particles. A rod-like particle interacts with
!! the wall through two embedded LJ sites.
!! 
!! @see Micheletti et. al. Journal of Chemical Physics 123, 224705.
!!
module m_lj2wall_interaction
  use m_particle, only: particle, single_interaction 
  use num_kind
  use class_poly_box
  use ljcylinder
  use json_module
  use m_json_wrapper, only: get_parameter
  use m_rod, only: rod
  use m_point, only: point
  use m_lj1wall_interaction, only: lj1wall_interaction
  implicit none

  type, extends(single_interaction) :: lj2wall_interaction
     type(lj1wall_interaction) :: a, b
     !> The distance of LJ interaction sites from the Gay-Berne particle
     !! center along the unique axis.
     real(dp) :: LJdist = 1.7
     
     !> If true, the wall interaction favors uniform alignment of rods 
     !! along the z-axis of the simulation box.
     logical :: isuniformalignment = .false.
     
   contains

     !> Computes the potential energy of the interaction between
     !! a particle and the wall.
     procedure :: potential => lj2wall_ia_potential
     
     !> Computes the force acting on the given particle by the wall.
     procedure :: force => lj2wall_ia_force

     !> Serializes the interaction to JSON.
     procedure :: to_json => lj2wall_ia_to_json

     !> Computes the interaction energy of a rod with the wall.
     procedure :: rod_potential => gbwall

     !> Returns the force acting on a rod by the wall.
     procedure :: rod_force => gbwall_force

     !> Returns the distances of the interactions sites in a rod
     !! from the cylinder axis.
     procedure :: rarb

     !> Returns the positions of the lj sites in a rod.
     procedure :: set_sites
 
     !> Produces a sample of the potential energies with this
     !! interaction.
     procedure :: sample => lj2wall_sample    
 end type lj2wall_interaction

 interface lj2wall_interaction
    module procedure lj2wall_interaction_from_json
 end interface lj2wall_interaction
 
contains

  
  !> Creates the lj2wall_interaction using parameters read from JSON
  !! contained in the @p json_val.
  function lj2wall_interaction_from_json(json_val) result(res)
    type(json_value), pointer, intent(in) :: json_val
    type(lj2wall_interaction) :: res
    call get_parameter(json_val, 'sig', res%a%sig, error_lb=0._dp)
    call get_parameter(json_val, 'wall_density', res%a%wall_density, &
         error_lb=0._dp)
    call get_parameter(json_val, 'eps', res%a%eps, error_lb=0._dp)
    res%b = res%a
    call get_parameter(json_val, 'alpha_A', res%a%alpha, error_lb=0._dp) 
    call get_parameter(json_val, 'alpha_B', res%b%alpha, error_lb=0._dp) 
    call get_parameter(json_val, 'LJ_dist', res%LJdist, error_lb=0._dp) 
    call get_parameter(json_val, 'is_uniform_alignment', res%isuniformalignment)
  end function


  !> Writes the parameters of this module with the format and output
  !! unit defined by @p writer.
  subroutine lj2wall_ia_to_json(this, json_val)
    class(lj2wall_interaction), intent(in) :: this
    type(json_value), intent(inout), pointer :: json_val
    
    call json_add(json_val, 'type', 'lj2wall_interaction')
    call json_add(json_val, 'wall_density', this%a%wall_density)
    call json_add(json_val, 'sig', this%a%sig)
    call json_add(json_val, 'eps', this%a%eps)

    call json_add(json_val, 'alpha_A', this%a%alpha) 
    call json_add(json_val, 'alpha_B', this%b%alpha) 
    call json_add(json_val, 'LJ_dist', this%LJdist) 
    call json_add(json_val, 'is_uniform_alignment', this%isuniformalignment)
  end subroutine


  !> Calculates the potential @p energy for a @p particle inside a
  !! cylindrical Lennard-Jones cavity. The radius of the cavity is
  !! defined by @p simbox. If the particle is too close to the wall or
  !! inside the wall @p overlap == true.
  pure subroutine lj2wall_ia_potential(this, particlei, simbox, energy, err)
    class(lj2wall_interaction), intent(in) :: this
    class(particle), intent(in) :: particlei
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    err = 0
    select type (particlei)
    type is (rod)
      call this%rod_potential(particlei, simbox, energy, overlap)
    class default
      err = 66
   end select
   if (overlap) err = 1
 end subroutine


  !> Returns the force acting on @p particlei by the wall of a
  !! cylindrical cavity. Radius of the cavity is defined by @p simbox.
  !! See module description for details.
  pure function lj2wall_ia_force(this, particlei, simbox) result(f)
    class(lj2wall_interaction), intent(in) :: this
    class(particle), intent(in) :: particlei
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)
    f = 0.
    select type (particlei)
    type is (rod)
      f = this%rod_force(particlei, simbox)
    end select
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
    class(lj2wall_interaction), intent(in) :: this
    class(rod), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    real(dp) :: rsite_a, rsite_b
    real(dp) :: fu, r_cylinder
    type(particle) :: site_a, site_b
    real(dp) :: energy_a, energy_b
    integer :: err_b, err_a
    overlap = .false.
    energy = 0._dp
    r_cylinder = getx(simbox)/2._dp 
    call this%set_sites(gbparticle, site_a, site_b)
    call this%a%potential(site_a, simbox, energy_a, err_a)
    call this%b%potential(site_b, simbox, energy_b, err_b)
    if (err_a == 0 .and. err_b == 0) then
       energy = energy_a + energy_b
       if (this%isuniformalignment) energy = angular(gbparticle) * energy
    else
       overlap = .true.
    end if
  end subroutine gbwall


  !> Returns the force exerted on a GB particle @p gbparticle by the
  !! wall of a cylindrical cavity. The @p gbparticle interacts with the
  !! wall via two embedded LJ sites.
  pure function gbwall_force(this, gbparticle, simbox) result(f)
    class(lj2wall_interaction), intent(in) :: this
    class(rod), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)
    real(dp) :: f_a(3), f_b(3)
    real(dp) :: rsite_a, rsite_b
    real(dp) :: fu, r_cylinder
    type(particle) :: site_a, site_b
    r_cylinder = getx(simbox) / 2. 
    call this%set_sites(gbparticle, site_a, site_b)
    f_a = this%a%force(site_a, simbox)
    f_b = this%b%force(site_b, simbox)
    f = f_a + f_b
    if (this%isuniformalignment) f = angular(gbparticle) * f
  end function


  !> Returns the distances from the cavity axis @p ra and @p ra for the
  !! interaction sites A and B, respectively. The interaction sites are
  !! embedded in @p arod.
  pure subroutine rarb(this, arod, ra, rb)
    class(lj2wall_interaction), intent(in) :: this
    class(rod), intent(in) :: arod
    real(dp), intent(out) :: ra, rb
    real(dp) :: xa, ya, xb, yb
    xa = arod%x + this%LJdist * arod%ux
    ya = arod%y + this%LJdist * arod%uy
    xb = arod%x - this%LJdist * arod%ux
    yb = arod%y - this%LJdist * arod%uy
    ra = sqrt(xa**2 + ya**2)
    rb = sqrt(xb**2 + yb**2)
  end subroutine


  !> Returns the distances from the cavity axis @p ra and @p ra for the
  !! interaction sites A and B, respectively. The interaction sites are
  !! embedded in @p arod.
  pure subroutine set_sites(this, arod, site_a, site_b)
    class(lj2wall_interaction), intent(in) :: this
    class(rod), intent(in) :: arod
    type(particle), intent(out) :: site_a, site_b
    site_a%x = arod%x + this%LJdist * arod%ux
    site_a%y = arod%y + this%LJdist * arod%uy
    site_b%x = arod%x - this%LJdist * arod%ux
    site_b%y = arod%y - this%LJdist * arod%uy    
  end subroutine set_sites

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
  subroutine lj2wall_sample(this, json_val, r, simbox)
    class(lj2wall_interaction), intent(in) :: this
    type(json_value), pointer :: json_val
    real(dp), intent(in) :: r(:)
    type(poly_box), intent(in) :: simbox
    type(rod) :: rodsamples(3)
    type(json_value), pointer :: child
    real(dp) :: energy
    integer :: err, i, k
    type(json_value), pointer :: r_json, energy_json
    character(len=3, kind=CK), parameter :: descriptions(3) = &
         ['GBx', 'GBy', 'GBz']

    ! First the GB-Wall interactions
    call rodsamples(1)%set_orientation([1._dp, 0._dp, 0._dp])
    call rodsamples(2)%set_orientation([0._dp, 1._dp, 0._dp])
    call rodsamples(3)%set_orientation([0._dp, 0._dp, 1._dp])
    call json_create_object(json_val, 'lj2wall_interaction')
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
    
  end subroutine lj2wall_sample

 
end module
