!> Contains the particlegroup type and related procedures.
module m_particlegroup
  use iso_fortran_env, only: dp => REAL64, output_unit, error_unit
  use class_poly_box, only: poly_box, minimage, isxperiodic, isyperiodic, &
       iszperiodic, getx, gety, getz
  use m_point, only: point, &
       movepoint_2, pair_interaction, pair_interaction_ptr, &
       single_interaction, single_interaction_ptr, &
       particlearray_wrapper, wrapper_delete, particlearray_to_json
  use utils, only: splitstr, join, acceptchange
  !$ use omp_lib
  use class_simplelist, only: simplelist, new_simplelist, simplelist_update, &
       simplelist_nbr_cells, flat_index, simplelist_deallocate, &
       simplelist_nbrmask, simplelist_cell_nbrmask
  use mt_stream, only: rng=>genrand_double1_s, rngstate=>mt_state
  use json_module
  use m_json_wrapper, only: get_parameter
  implicit none  

  !> The particlegroup to simplify handling multiple cell lists and and
  !! arrays of particles.
  type particlegroup
     !> The name of the group
     character(len=:), allocatable :: name
     !> The particles in this group.
     class(point), allocatable :: particles(:)
     !> The cell list
     type(simplelist) :: sl
   contains
     procedure :: to_json => particlegroup_to_json
     procedure :: scalepositions
     procedure :: check_particles => particlegroup_check_particles
     final :: particlegroup_finalize
  end type particlegroup

  !> Wrapper to be used for arrays of particlegroups of different types.
  type particlegroup_ptr
     class(particlegroup), pointer :: ptr => null()
  end type particlegroup_ptr

  !> Constructor interface.
  interface particlegroup
     module procedure create_particlegroup
  end interface particlegroup
  
contains

  
  !> Creates a particlegroup.
  !!
  !! @param simbox the box in which the particles reside.
  !! @param particles the particles in the group.
  !! @param min_cell_length the minimum cell side length in the cell
  !!        list.
  !! @param min_boundary_width the width of the "boundary" area around
  !!        the cell. When any particle crosses the boundary, the cell
  !!        list should be updated.
  !! 
  !! @return the particlegroup which has a cell list with the given
  !!         specs.
  !!
  function create_particlegroup(simbox, particles, min_cell_length, &
       min_boundary_width, name) result(group)
    type(poly_box), intent(in) :: simbox
    class(point), intent(in) :: particles(:)
    real(dp), intent(in) :: min_cell_length, min_boundary_width
    character(kind=CK, len=*), intent(in) :: name
    type(particlegroup) :: group
    group%name = name
    allocate(group%particles(size(particles)), &
         source=particles(1:size(particles)))
    !$ if (.true.) then
    !$ call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
    !$& min_boundary_width, is_nx_even = isxperiodic(simbox), &
    !$& is_ny_even = isyperiodic(simbox), is_nz_even = iszperiodic(simbox))
    !$ else 
    call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
         min_boundary_width)
    !$ end if
  end function create_particlegroup


  !> Serializes the particlegroup @p this to JSON @p json_val.
  subroutine particlegroup_to_json(this, json_val)
    class(particlegroup), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    if (allocated(this%name)) call json_add(json_val, 'name', this%name)
    !! particlearray
    call particlearray_to_json(json_val, this%particles)
  end subroutine particlegroup_to_json

  
  !> Final routine.
  subroutine particlegroup_finalize(group)
    type(particlegroup), intent(inout) :: group
    call simplelist_deallocate(group%sl)
    if(allocated(group%particles)) deallocate(group%particles)
  end subroutine particlegroup_finalize
  

  !> Computes the total @p energy of the system. Interactions are
  !! computed cell-by-cell and group-by-group.
  !!
  !! @param groups the particlegroups to which the total_energy is
  !!        computed.
  !! @param simbox the box in which the groups reside.
  !! @param pair_ias the matrix of pair_interactions.
  !! @param single_ias the array of single_interactions.
  !! @param energy the total energy at return.
  !! @param err is non-zero if an error occurs, such as two particles
  !!        being too close to each other.
  !!
  subroutine total_energy(groups, simbox, pair_ias, single_ias, energy, &
       err)
    type(particlegroup_ptr), intent(inout) :: groups(:)
    type(poly_box), intent(in) :: simbox
    type(pair_interaction_ptr), intent(in) :: pair_ias(:, :)
    type(single_interaction_ptr), intent(in) :: single_ias(:)
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    integer :: i, ix, iy, iz, i_group
    real(dp) :: energy_j
    integer :: nbr_cells(3, 27), n_nbr_cells
    integer :: j_group
    energy = 0._dp
    err = 0
    if (size(groups) == 0) return
    !$OMP PARALLEL default(shared) reduction(+:energy, err)& 
    !$OMP& private(energy_j, i, nbr_cells, n_nbr_cells)
    !$OMP DO collapse(3) schedule(dynamic)
    do ix = 0, groups(1)%ptr%sl%nx - 1 
       do iy = 0, groups(1)%ptr%sl%ny - 1
          do iz = 0, groups(1)%ptr%sl%nz - 1
             if (err == 0) then
                do i_group = 1, size(groups)
                   ! 1. compute inside ix, iy, iz in i_group
                   call cell_energy(groups(i_group)%ptr, ix, iy, iz, simbox, &
                        pair_ias(i_group, i_group)%ptr, &
                        single_ias(i_group)%ptr, energy_j, err)
                   if (err /= 0) exit
                   energy = energy + energy_j
                end do
             end if
             if (err == 0) then
                do i_group = 1, size(groups) - 1
                   !! 2. compute with ix, iy, iz in j_group > i_group
                   do j_group = i_group + 1, size(groups)
                      !write(*, *) "compute with ", ix, iy, iz, " in ",
                      ! j_group, " > ", i_group
                      call cell_pair_energy(groups(i_group)%ptr, ix, iy, iz, &
                           groups(j_group)%ptr, ix, iy, iz, simbox, &
                           pair_ias(i_group, j_group)%ptr, energy_j, err)
                      if (err /= 0) exit
                      energy = energy + energy_j
                   end do
                   if (err /= 0) exit
                end do
             end if
             if (err == 0) then
                ! 3. for all j_group (including i_group) compute where
                ! ix + nx * iy + nx * ny * iz < jx + jy * nx + jz * nx * ny
                ! and jx, jy, jz is a neighbour of ix, iy, iz.
                call simplelist_nbr_cells(groups(i_group)%ptr%sl, ix, iy, iz, &
                     nbr_cells, n_nbr_cells)
                do i_group = 1, size(groups)
                   do i = 1, n_nbr_cells
                      if (flat_index(groups(i_group)%ptr%sl, nbr_cells(1, i), &
                           nbr_cells(2, i), nbr_cells(3, i)) > &
                           flat_index(groups(i_group)%ptr%sl, ix, iy, iz)) then
                         do j_group = 1, size(groups)
                            call cell_pair_energy(groups(i_group)%ptr, &
                                 ix, iy, iz, &
                                 groups(j_group)%ptr, nbr_cells(1, i), &
                                 nbr_cells(2, i), nbr_cells(3, i), simbox, &
                                 pair_ias(i_group, j_group)%ptr, energy_j, err)
                            if (err /= 0) exit
                            energy = energy + energy_j
                         end do
                      end if
                      if (err /= 0) exit
                   end do
                   if (err /= 0) exit
                end do
             end if
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL  
  end subroutine total_energy

  !> The "internal" energy of cell @p ix, @p iy, @p iz in @p this group.
  !!
  !! @param this the particlegroup
  !! @param ix, iy, iz the cell indices.
  !! @param simbox the simulation box.
  !! @param pair_ia the pair interaction between particles in this
  !!        group.
  !! @param single_ia the single interaction concerning this group.
  !! @param energy the energy at return.
  !! @param err is non-zero if any error occurs.
  !!
  subroutine cell_energy(this, ix, iy, iz, simbox, pair_ia, single_ia, &
       energy, err)
    type(particlegroup), intent(in) :: this
    integer, intent(in) :: ix, iy, iz
    type(poly_box), intent(in) :: simbox
    class(pair_interaction), intent(in) :: pair_ia
    class(single_interaction), pointer, intent(in) :: single_ia
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    integer :: i, j
    real(dp) :: energy_ij, rij(3)
    energy = 0
    err = 0
    do i = 1, this%sl%counts(ix, iy, iz) - 1
       associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
         do j = i + 1, this%sl%counts(ix, iy, iz)
            associate(particlej => &
                 this%particles(this%sl%indices(j, ix, iy, iz)))
              rij = minimage(simbox,&
                   particlej%x - particlei%x, particlej%y - particlei%y,&
                   particlej%z - particlei%z)
              if (dot_product(rij, rij) < pair_ia%get_cutoff()**2) then
                 call pair_ia%pair_potential(particlei, particlej, rij, &
                      energy_ij, err)
                 if (err /= 0 ) return
                 energy = energy + energy_ij
              end if
            end associate
         end do
       end associate
    end do
    if (associated(single_ia)) then
       do i = 1, this%sl%counts(ix, iy, iz)
          associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
            call single_ia%potential(particlei, simbox, energy_ij, err)
            if (err /= 0) exit
            energy = energy + energy_ij
          end associate
       end do
    end if
  end subroutine cell_energy

  
  !> The total pair interaction between particles in groups @p this and
  !! @p another in cells [ix, iy, iz] and [jx, jy, jz], respectively.
  !!
  !! @param this the particlegroup
  !! @param ix, iy, iz the cell index for @p this.
  !! @param another the other particlegroup.
  !! @param jx, jy, jz the cell index for @p another.
  !! @param simbox the simulation box in which the particles reside.
  !! @param pair_ia the pair_interaction between @p this and @p another.
  !! @param energy the energy at return.
  !! @param err is non-zero if any pairwise-computed interaction
  !!        results in non-zero err.
  !!
  subroutine cell_pair_energy(this, ix, iy, iz, another, jx, jy, jz, simbox, &
     pair_ia, energy, err)
  type(particlegroup), intent(in) :: this, another
  integer, intent(in) :: ix, iy, iz, jx, jy, jz
  type(poly_box), intent(in) :: simbox
  class(pair_interaction), intent(in) :: pair_ia
  real(dp), intent(out) :: energy
  integer, intent(out) :: err
  integer :: i, j
  real(dp) :: energy_ij, rij(3)
  !write(*, *) "compute ", this%name, ix, iy, iz, " with ", &
  !     another%name, jx, jy, jz
  energy = 0
  if ((ix-jx)**2 + (iy-jy)**2 + (iz-jz)**2 <= 3) then
  ! No need to compute minimum image.
  do i = 1, this%sl%counts(ix, iy, iz)
     associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
       do j = 1, another%sl%counts(jx, jy, jz)
          associate(particlej => another%particles(&
               another%sl%indices(j, jx, jy, jz)))
            rij = [particlej%x - particlei%x, particlej%y - particlei%y, &
                 particlej%z - particlei%z]
            if (dot_product(rij, rij) < pair_ia%get_cutoff()**2) then
               call pair_ia%pair_potential(particlei, particlej, rij, &
                    energy_ij, err)
               if (err /= 0) return
               energy = energy + energy_ij
            end if
          end associate
       end do
     end associate
  end do
  else
  do i = 1, this%sl%counts(ix, iy, iz)
     associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
       do j = 1, another%sl%counts(jx, jy, jz)
          associate(particlej => another%particles(&
               another%sl%indices(j, jx, jy, jz)))
            rij = minimage(simbox, particlej%x - particlei%x, &
                 particlej%y - particlei%y, particlej%z - particlei%z)
            if (dot_product(rij, rij) < pair_ia%get_cutoff()**2) then
               call pair_ia%pair_potential(particlei, particlej, rij, &
                    energy_ij, err)
               if (err /= 0) return
               energy = energy + energy_ij
            end if
          end associate
       end do
     end associate
  end do
  end if
end subroutine cell_pair_energy


!> Scales the positions of @p particles with the same factors that were
!! used to scale the simulation box dimensions from @p oldbox to
!! @p newbox.
!!
!! @param this the particlegroup for which the particles are scaled.
!! @param oldbox, newbox the new and old simulation boxes.
!! 
subroutine scalepositions(this, oldbox, newbox)
  class(particlegroup), intent(inout) :: this
  type(poly_box), intent(in) :: oldbox
  type(poly_box), intent(in) :: newbox
  integer :: i
  !! Don't turn the loop below to array syntax. At least GCC6.0.0 (BETA) can
  !! not handle these arrays correctly.
  do i = 1, size(this%particles)
     this%particles(i)%x = this%particles(i)%x * getx(newbox) / getx(oldbox)
     this%particles(i)%y = this%particles(i)%y * gety(newbox) / gety(oldbox)
     this%particles(i)%z = this%particles(i)%z * getz(newbox) / getz(oldbox)
  end do
  call simplelist_update(this%sl, newbox, this%particles)
end subroutine scalepositions

!> Checks that all particles in @p this group are inside the @p simbox.
!!
!! @param the particlegroup to check.
!! @param simbox the simulation box in which the particles should
!!        reside.
function particlegroup_check_particles(this, simbox) result(res)
  class(particlegroup), intent(in) :: this
  class(poly_box), intent(in) :: simbox
  logical :: res
  integer :: i
  res = .true.
  do i = 1, size(this%particles)
     res = res .and. simbox%check_position(this%particles(i)%position())
  end do
end function particlegroup_check_particles

end module m_particlegroup
