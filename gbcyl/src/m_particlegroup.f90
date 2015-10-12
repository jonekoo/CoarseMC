module m_particlegroup
  use particle, only: particledat
  use class_simplelist, only: simplelist, new_simplelist, simplelist_deallocate
  use energy, only: get_cutoff
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

  subroutine init(group, simbox)
    class(particlegroup) :: group
    type(poly_box), intent(in) :: simbox
    real(dp) :: min_cell_length
    min_cell_length = get_cutoff() + 2._dp * get_max_translation()
    !$ if (.true.) then
    !$ call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
    !$& is_x_even = isxperiodic(simbox), is_y_even = isyperiodic(simbox), &
    !$& is_z_even = iszperiodic(simbox), cutoff=get_cutoff())
    !$ else 
    call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
         cutoff=get_cutoff())
    !$ end if
  end subroutine init

  impure elemental subroutine particlegroup_finalize(group)
    type(particlegroup), intent(inout) :: group
    call simplelist_deallocate(group%sl)
    if(allocated(group%particles)) deallocate(group%particles)
  end subroutine particlegroup_finalize
  
end module m_particlegroup
