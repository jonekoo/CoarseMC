!> Implements Metropolis Monte Carlo updates of the coordinates of the
!! particles and the simulation box dimensions. The module also
!! controls the parallel tempering updates. If compiled and run
!! with OpenMP, the domain decomposition algorithm is used for the
!! particle moves and energy calculations. 
module m_particlegroup
  use iso_fortran_env, only: dp => REAL64, output_unit, error_unit
  use class_poly_box, only: poly_box, minimage, isxperiodic, isyperiodic, &
       iszperiodic
  use particle, only: particledat, position, setposition, &
       moveparticle_2, pair_interaction, pair_interaction_ptr, &
       single_interaction, single_interaction_ptr, &
       particlearray_wrapper, wrapper_delete, particlearray_to_json
  use class_parameterizer, only: parameterizer, getparameter
  use class_parameter_writer, only: parameter_writer, writeparameter, &
       writecomment
  use genvoltrial
  use utils, only: splitstr, join, acceptchange
  !$ use omp_lib
  use class_simplelist, only: simplelist, new_simplelist, simplelist_update, &
       simplelist_nbr_cells, flat_index, simplelist_deallocate, &
       simplelist_nbrmask, simplelist_cell_nbrmask
  include 'rng.inc'
  use json_module
  use m_json_wrapper, only: get_parameter
  implicit none  

  type particlegroup
     character(len=:), allocatable :: name
     class(particledat), allocatable :: particles(:)
     type(simplelist) :: sl
   contains
     procedure :: to_json => particlegroup_to_json
     procedure :: scalepositions
     final :: particlegroup_finalize
  end type particlegroup

  type particlegroup_ptr
     class(particlegroup), pointer :: ptr => null()
  end type particlegroup_ptr

  interface particlegroup
     module procedure create_particlegroup
  end interface particlegroup

contains

  function create_particlegroup(simbox, particles, min_cell_length, &
       min_boundary_width, name) result(group)
    type(poly_box), intent(in) :: simbox
    class(particledat), intent(in) :: particles(:)
    real(dp), intent(in) :: min_cell_length, min_boundary_width
    character(kind=CK, len=*), intent(in) :: name
    type(particlegroup) :: group
    group%name = name
    allocate(group%particles(size(particles)), source=particles)
    !$ if (.true.) then
    !$ call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
    !$& min_boundary_width, is_nx_even = isxperiodic(simbox), &
    !$& is_ny_even = isyperiodic(simbox), is_nz_even = iszperiodic(simbox))
    !$ else 
    call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
         min_boundary_width)
    !$ end if
  end function create_particlegroup

  subroutine particlegroup_to_json(this, json_val)
    class(particlegroup), intent(in) :: this
    type(json_value), pointer, intent(inout) :: json_val
    if (allocated(this%name)) call json_add(json_val, 'name', this%name)
    !! particlearray
    call particlearray_to_json(json_val, this%particles)
    !! cell list

  end subroutine particlegroup_to_json
  
  impure elemental subroutine particlegroup_finalize(group)
    type(particlegroup), intent(inout) :: group
    call simplelist_deallocate(group%sl)
    if(allocated(group%particles)) deallocate(group%particles)
  end subroutine particlegroup_finalize
  

!> Computes the total @p energy of the system. Interactions are computed
!! cell-by-cell for the @p particles in the cell list @p sl. @p simbox
!! is the simulation cell. If any two particles are too close to each
!! other, @p overlap is true. 
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
                 !! 1. compute inside ix, iy, iz in i_group
                 call cell_energy(groups(i_group)%ptr, ix, iy, iz, simbox, &
                      pair_ias(i_group, i_group)%ptr, single_ias(i_group)%ptr, &
                      energy_j, err)
                 if (err /= 0) exit
                 energy = energy + energy_j
              end do
           end if
           if (err == 0) then
              do i_group = 1, size(groups) - 1
                 !! 2. compute with ix, ix, y in j_group > i_group
                 do j_group = i_group + 1, size(groups)
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
              !! 3. for all j_group (including i_group) compute where
              !! ix + nx * iy + nx * ny * iz < jx + jy * nx + jz * nx * ny
              !! and jx, jy, jz is a neighbour of ix, iy, iz.
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
        associate(particlej => this%particles(this%sl%indices(j, ix, iy, iz)))
          rij = minimage(simbox,&
               particlej%x - particlei%x, particlej%y - particlei%y,&
               particlej%z - particlei%z)
          if (norm2(rij) < pair_ia%get_cutoff()) then
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

  energy = 0
  do i = 1, this%sl%counts(ix, iy, iz)
     associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
       do j = 1, another%sl%counts(jx, jy, jz)
          associate(particlej => another%particles(&
               another%sl%indices(j, jx, jy, jz)))
            rij = minimage(simbox, particlej%x - particlei%x, &
                 particlej%y - particlei%y, particlej%z - particlei%z)
            if (norm2(rij) < pair_ia%get_cutoff()) then
               call pair_ia%pair_potential(particlei, particlej, rij, &
                    energy_ij, err)
               if (err /= 0) return
               energy = energy + energy_ij
            end if
          end associate
     end do
     end associate
  end do
end subroutine cell_pair_energy


!> Scales the positions of @p particles with the same factors that were
!! used to scale the simulation box dimensions from @p oldbox to
!! @p newbox.
subroutine scalepositions(this, oldbox, newbox)
  class(particlegroup), intent(inout) :: this
  type(poly_box), intent(in) :: oldbox
  type(poly_box), intent(in) :: newbox
  this%particles(:)%x = this%particles(:)%x * getx(newbox) / getx(oldbox)
  this%particles(:)%y = this%particles(:)%y * gety(newbox) / gety(oldbox)
  this%particles(:)%z = this%particles(:)%z * getz(newbox) / getz(oldbox)
  call simplelist_update(this%sl, newbox, this%particles)
end subroutine scalepositions

end module m_particlegroup
