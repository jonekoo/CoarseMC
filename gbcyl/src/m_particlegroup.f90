module m_particlegroup
  use particle, only: particledat
  use class_simplelist, only: simplelist, new_simplelist, simplelist_deallocate
  use particle_mover, only: get_max_translation
  use class_parameterizer, only: parameterizer
  use class_poly_box
  use nrtype, only: dp
  implicit none
  
  type particlegroup
     type(particledat), allocatable :: particles(:)
     type(simplelist) :: sl
   contains
     procedure :: init
     final :: particlegroup_finalize
  end type particlegroup

contains

  subroutine init(group, simbox, min_cell_length, min_boundary_width)
    class(particlegroup) :: group
    type(poly_box), intent(in) :: simbox
    real(dp), intent(in) :: min_cell_length, min_boundary_width
    !$ if (.true.) then
    !$ call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
    !$& min_boundary_width, is_x_even = isxperiodic(simbox), &
    !$& is_y_even = isyperiodic(simbox), is_z_even = iszperiodic(simbox))
    !$ else 
    call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
         min_boundary_width)
    !$ end if
  end subroutine init

  impure elemental subroutine particlegroup_finalize(group)
    type(particlegroup), intent(inout) :: group
    call simplelist_deallocate(group%sl)
    if(allocated(group%particles)) deallocate(group%particles)
  end subroutine particlegroup_finalize
  
end module m_particlegroup
