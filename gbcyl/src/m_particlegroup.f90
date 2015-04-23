module m_particlegroup
  use m_particle
  use class_poly_box
  use m_simplelist
  use particle_mover
  use energy, only: get_cutoff
  implicit none
  
  type particlegroup
  end type particlegroup

  type, extends(particlegroup) :: polymorphic_group
     class(particle), allocatable :: particles(:)
     type(simplelist) :: sl
     type(poly_box), pointer :: simbox => null()
   contains
     procedure :: update => polymorphic_update
  end type polymorphic_group

  interface polymorphic_group
     procedure polymorphic_from_particles
  end interface polymorphic_group
  
contains

  
  
  function polymorphic_from_particles(particles, simbox) result(this)
    class(particle), intent(in) :: particles(:)
    type(poly_box), pointer, intent(in) :: simbox
    type(polymorphic_group) :: this
    real(dp) :: min_cell_length
    this%simbox => simbox
    allocate(this%particles(size(particles)), source=particles)
    min_cell_length = get_cutoff() + 2 * get_max_translation()
    !$ if (.true.) then
    !$ call new_simplelist(this%sl, this%simbox, this%particles, &
    !$& min_cell_length, &
    !$& is_x_even = isxperiodic(simbox), is_y_even = isyperiodic(simbox), &
    !$& is_z_even = iszperiodic(simbox), cutoff=get_cutoff())
    !$ else 
    call new_simplelist(this%sl, this%simbox, this%particles, min_cell_length, &
         cutoff=get_cutoff())
    !$ end if
  end function polymorphic_from_particles

  subroutine polymorphic_update(this)
    class(polymorphic_group), intent(inout) :: this
    call update(this%sl, this%simbox, this%particles)
  end subroutine polymorphic_update
    
  
end module m_particlegroup
