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

  !> The maximum absolute change of volume in a trial volume update.
  real(dp), save :: maxscaling = 100._dp

  !> The types of trial volume updates.
  character(len = 200), dimension(:), allocatable, save :: scalingtypes

  !> True if the module is correctly initialized.
  logical :: is_initialized = .false.

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

  type, extends(particlearray_wrapper) :: domain
    integer :: n_cell = 0
  contains
    procedure :: domain_assign
    generic :: assignment(=) => domain_assign
    final :: domain_delete
  end type

  interface particlegroup_init
     module procedure mcsweep_from_json, mcsweep_init
  end interface particlegroup_init
  
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
  
  impure elemental subroutine domain_assign(this, src)
    class(domain), intent(inout) :: this
    type(domain), intent(in) :: src
    this%particlearray_wrapper = src%particlearray_wrapper
    this%n_cell = src%n_cell
  end subroutine domain_assign

  impure elemental subroutine domain_delete(this)
    type(domain), intent(inout) :: this
    this%n_cell = 0
  end subroutine domain_delete
  
!> Initializes the module by getting the module parameters from the
!! @p reader object.
!!
!! @param reader the object which gets parameters by name e.g. from an
!!        input file.
!! @param the_simbox the simulation box defining the boundary
!!        conditions and volume.
!! @param the_particles the particles to be updated with the MC moves.
!!
subroutine mcsweep_init(reader)
  type(parameterizer), intent(in) :: reader
  character(len = 200), save :: scalingtype = "z"
  call getparameter(reader, 'scaling_type', scalingtype)
  call parsescalingtype(scalingtype)
  call getparameter(reader, 'max_scaling', maxscaling)
  is_initialized = .true.
end subroutine mcsweep_init

subroutine mcsweep_from_json(json_val)
  type(json_value), pointer, intent(in) :: json_val
  allocate(scalingtypes(0))
  call get_parameter(json_val, 'scaling_types', scalingtypes) 
  call get_parameter(json_val, 'max_scaling', maxscaling, error_lb=0._dp)
  is_initialized = .true.
end subroutine mcsweep_from_json

!> Finalizes the module state.
subroutine mcsweep_finalize
  if (allocated(scalingtypes)) deallocate(scalingtypes)
  is_initialized = .false.
end subroutine mcsweep_finalize

!> Checks and parses the @p scalingtype string to the scalingtypes
!! array. 
subroutine parsescalingtype(scalingtype)
  character(len=*), intent(in) :: scalingtype
  logical :: isok
  isok = (0 == verify(trim(adjustl(scalingtype)), 'xyz,') .and. &
       0 /= len_trim(adjustl(scalingtype)))
  if (.not. isok) then
     write(error_unit, *) 'Error parsing parameter scaling_type, stopping!'
     stop 'Stopped by parsescalingtype'
  end if
  call splitstr(scalingtype, ',', scalingtypes)
end subroutine parsescalingtype

subroutine set_maxscaling(upper_limit)
  real(dp), intent(in) :: upper_limit
  maxscaling = upper_limit
end subroutine set_maxscaling

function get_maxscaling() result(upper_limit)
  real(dp) :: upper_limit
  upper_limit = maxscaling
end function get_maxscaling

!> Writes the parameters and observables of this module and its
!! dependencies. The @p writer defines the format and output unit..
subroutine mc_sweep_writeparameters(writer)
  type(parameter_writer), intent(inout) :: writer
  character(len=200) :: joined
  call writecomment(writer, 'MC sweep parameters')
  call writeparameter(writer, 'max_scaling', maxscaling)
  if (allocated(scalingtypes)) then
     call join(scalingtypes, ',', joined)
  else
     joined = ''
  end if
  call writeparameter(writer, 'scaling_type', trim(joined))
end subroutine mc_sweep_writeparameters

!> Writes the parameters and observables of this module and its
!! dependencies. The @p writer defines the format and output unit..
subroutine mc_sweep_to_json(json_val)
  type(json_value), intent(inout), pointer :: json_val
  type(json_value), pointer :: temp
  type(json_value), pointer :: str
  integer :: i
  if (allocated(scalingtypes)) then
     call json_add(json_val, 'max_scaling', maxscaling)
     if (size(scalingtypes) > 0) then
        call json_create_array(temp, 'scaling_types')
        do i = 1, size(scalingtypes)
           call json_create_string(str, scalingtypes(i), '')
           call json_add(temp, str)
        end do
        call json_add(json_val, temp)
     end if
  else
     stop 'ERROR: scalingtypes not allocated!'
  end if
end subroutine mc_sweep_to_json

!> Schedules parallel moves of particles using OpenMP with a domain 
!! decomposition algorithm. 
!!
!! @see e.g. G. Heffelfinger and M. Lewitt. J. Comp. Chem., 17(2):250–265,
!! 1996.
!! 
subroutine make_particle_moves(groups, genstates, simbox, temperature, &
     pair_ias, single_ias, dE, n_trials, n_accepted)
  type(particlegroup_ptr), intent(inout) :: groups(:)
  type(poly_box), intent(in) :: simbox
  type(rngstate), intent(inout) :: genstates(0:)
  type(pair_interaction_ptr), intent(in) :: pair_ias(:, :)
  type(single_interaction_ptr) :: single_ias(:)
  real(dp), intent(in) :: temperature
  real(dp), intent(out) :: dE
  integer, intent(out) :: n_accepted, n_trials
  !$ integer :: n_threads
  integer :: thread_id
  !! You may be tempted to make the allocatable arrays automatic, but there's
  !! no performance gained and depending on the system large automatic arrays
  !! may cause a stack overflow.
  real(dp) :: dE_d
  integer :: ix, iy, iz, jx, jy, jz, n_trials_d, n_accepted_d, i_group, i
  type(domain), allocatable :: ds(:)
  do i_group = 1, size(groups)
     call simplelist_update(groups(i_group)%ptr%sl, simbox, &
          groups(i_group)%ptr%particles)
  end do
  thread_id = 0
  dE = 0._dp
  n_accepted = 0
  n_trials = 0
  if (size(groups) == 0) return
  !$ n_threads = 1
  !! Loop over cells. This can be thought of as looping through a 
  !! 2 x 2 x 2 cube of cells.
  !$OMP PARALLEL shared(groups, simbox, genstates, pair_ias, single_ias)& 
  !$OMP& private(thread_id, n_threads, n_trials_d, n_accepted_d, ds, dE_d, &
  !$OMP& i_group)&
  !$OMP& reduction(+:dE, n_accepted, n_trials) 
  !$ thread_id = omp_get_thread_num()
  allocate(ds(size(groups)))
  do iz=0, min(1, groups(1)%ptr%sl%nz-1)
     do iy=0, min(1, groups(1)%ptr%sl%ny-1)
        do ix=0, min(1, groups(1)%ptr%sl%nx-1)
           !$OMP DO collapse(3) schedule(dynamic)
           do jz = iz, groups(1)%ptr%sl%nz - 1, 2
              do jy = iy, groups(1)%ptr%sl%ny - 1, 2
                 do jx = ix, groups(1)%ptr%sl%nx - 1, 2

                    !! Collect temp_particles for all groups
                    do i_group = 1, size(groups)
                       ds(i_group) = create_domain(groups(i_group)%ptr, simbox,&
                            jx, jy, jz)
                    end do
                    !! Move particles
                    do i_group = 1, size(groups)
                       call domain_move(ds, i_group, &
                            genstates(thread_id:thread_id), simbox, &
                            temperature, pair_ias, &
                            single_ias, dE_d, n_trials=n_trials_d, &
                            n_accepted=n_accepted_d)
                       dE = dE + dE_d
                       n_trials = n_trials + n_trials_d
                       n_accepted = n_accepted + n_accepted_d
                    end do
                    !! Synchronize
                    !! :TODO: Generalize to a sync domain operation?
                    do i_group = 1, size(groups)
                       !! Cray compiler does not accept this since "An actual
                       !! argument must be definable when associated with a
                       !! dummy argument that has intent(out) or intent(inout)."
                       do i = 1, ds(i_group)%n_cell
                          call groups(i_group)%ptr%particles(&
                               groups(i_group)%ptr%sl%indices(i, jx, jy, jz)&
                               )%downcast_assign(ds(i_group)%arr(i))
                       end do
                    end do
                    call domain_delete(ds)
                 end do
              end do
           end do
           !$OMP END DO 
           !! The end of parallelized loop forces an implicit barrier.
           !! Memory view of the threads is also synchronized here.
        end do
        !$OMP BARRIER
     end do
     !$OMP BARRIER
  end do
  !$OMP END PARALLEL
  do i_group = 1, size(groups)
     call simplelist_update(groups(i_group)%ptr%sl, simbox, &
          groups(i_group)%ptr%particles)
  end do
end subroutine make_particle_moves


function create_domain(this, simbox, jx, jy, jz) result(d)
  type(particlegroup), intent(in) :: this
  type(poly_box), intent(in) :: simbox
  integer, intent(in) :: jx, jy, jz
  type(domain) :: d
  logical, allocatable :: nbr_mask(:)
  integer :: n_local
  integer :: i
  integer, allocatable :: helper(:)

  allocate(nbr_mask(size(this%particles)), source=.false.)
  d%n_cell = this%sl%counts(jx, jy, jz)
  call simplelist_cell_nbrmask(this%sl, simbox, jx, jy, jz, nbr_mask)
  n_local = count(nbr_mask)
  !! Remove particles in jx, jy, jz from mask:
  nbr_mask(this%sl%indices(1:d%n_cell, jx, jy, jz)) = .false.

  !! GFortran 6.0.0 is not fine with pack from a polymorphic
  !! array class(particledat) so we'll use a workaround.
  !! Put particles in jx, jy, jz first in temp_particles:
  allocate(helper, source=[this%sl%indices(1:d%n_cell, jx, jy, jz), &
       pack([(i, i = 1, size(this%particles))], nbr_mask)])

  allocate(d%arr(n_local), mold=this%particles(1))
  do i = 1, n_local
     call d%arr(i)%downcast_assign(this%particles(helper(i)))
  end do

  allocate(d%mask(n_local), source=.true.)
end function create_domain

impure elemental subroutine delete_domain(this)
  class(domain), intent(inout) :: this
  this%n_cell = 0
end subroutine delete_domain

subroutine domain_move(domains, i_d, genstates, simbox, temperature, &
     pair_ias, single_ias, dE, n_trials, n_accepted)
  type(domain), intent(inout) :: domains(:)
  integer, intent(in) :: i_d
  type(rngstate), intent(inout) :: genstates(:)
  type(poly_box), intent(in) :: simbox
  real(dp), intent(in) :: temperature
  type(pair_interaction_ptr), intent(in) :: pair_ias(:, :)
  type(single_interaction_ptr) :: single_ias(:)
  real(dp), intent(out) :: dE
  integer, intent(out), optional :: n_trials, n_accepted
  integer :: j
  class(particledat), allocatable :: newparticle
  integer :: err
  real(dp) :: enew
  real(dp) :: eold
  logical :: isaccepted
  n_accepted = 0
  n_trials = 0
  dE = 0.
  if (domains(i_d)%n_cell > 0) allocate(newparticle, source=domains(i_d)%arr(1))
  do j = 1, domains(i_d)%n_cell
     domains(i_d)%mask(j) = .false.
     call newparticle%downcast_assign(domains(i_d)%arr(j))
     call newparticle%move(genstates(1))
     call setposition(newparticle, minimage(simbox, newparticle%x, &
          newparticle%y, newparticle%z))
     call newparticle%energy(domains, pair_ias(:, i_d), simbox, &
          single_ias(i_d)%ptr, enew, err)
     if(err == 0) then 
        call domains(i_d)%arr(j)%energy(domains, pair_ias(:, i_d), simbox, &
             single_ias(i_d)%ptr, eold, err)
        if (err /= 0) then
           call domains(i_d)%arr(j)%to_stdout()
           call newparticle%to_stdout()
           write(error_unit, *) 'ERROR: err=', err, ' domain=', i_d, ' old particle ', j
           stop 'Stopped by domain_move.'
        end if
        call acceptchange(eold, enew, temperature, genstates(1), isaccepted)
        if(isaccepted) then
           call domains(i_d)%arr(j)%downcast_assign(newparticle)
           dE = dE + enew - eold
           n_accepted = n_accepted + 1
        end if
     end if
     domains(i_d)%mask(j) = .true.
  end do
  n_trials = n_trials + domains(i_d)%n_cell
end subroutine domain_move



subroutine update_volume(simbox, groups, genstate, pair_ias, &
     single_ias, temperature, pressure, e_total, n_trials, n_accepted)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup_ptr), intent(inout) :: groups(:)
  type(rngstate), intent(inout) :: genstate
  type(pair_interaction_ptr), intent(in) :: pair_ias(:, :)
  type(single_interaction_ptr) :: single_ias(:)
  real(dp), intent(in) :: temperature, pressure
  real(dp), intent(out) :: e_total
  integer, intent(out), optional :: n_trials, n_accepted
  integer :: i
  do i = 1, size(scalingtypes)
     call movevol(simbox, groups, scalingtypes(i), genstate, pair_ias, &
          single_ias, temperature, pressure, e_total, n_trials, &
          n_accepted)
  end do
end subroutine update_volume


!> Performs a trial volume scaling which scales the @p simbox and all
!! the positions of @p particles. For more information see for example
!! Allen and Tildesley: Computer Simulation of Liquids: the chapter
!! about NPT ensemble Monte Carlo.
!! 
!! @param simbox the simulation box.
!! @param particles the particles in the simulation box.
!! @param scalingtype string defining the direction(s) for changing
!!        the system volume.
!! @param genstate the random number generator state.
!! 
subroutine movevol(simbox, groups, scalingtype, genstate, pair_ias, &
     single_ias, temperature, pressure, e_total, n_trials, n_accepted)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup_ptr), intent(inout) :: groups(:)
  character(len=*), intent(in) :: scalingtype
  type(rngstate), intent(inout) :: genstate
  type(pair_interaction_ptr), intent(in) :: pair_ias(:, :)
  type(single_interaction_ptr) :: single_ias(:)
  real(dp), intent(in) :: temperature, pressure
  real(dp), intent(out) :: e_total
  integer, intent(out), optional :: n_trials, n_accepted
  integer :: nparticles, err
  real(dp) :: Vo, Vn
  logical :: isaccepted
  real(dp) :: totenew
  real(dp) :: boltzmannn
  real(dp) :: boltzmanno
  type(poly_box) :: oldbox
  real(dp), dimension(3) :: scaling
  integer :: i
  nparticles = 0
  do i = 1, size(groups)
     nparticles = nparticles + size(groups(i)%ptr%particles)
  end do
  !! Store old volume and simulation box
  Vo = volume(simbox)
  oldbox = simbox
  
  !! It seems that total energy may drift (WHY?!) if it is not updated here:
  call total_energy(groups, simbox, pair_ias, single_ias, &
       e_total, err)
  if (err /= 0) stop 'movevol: overlap in old configuration! '//&
       'Should never happen!'
  
  !! Scale coordinates and the simulations box
  scaling = genvoltrial_scale(simbox, maxscaling, genstate, &
       trim(adjustl(scalingtype)))
  do i = 1, size(groups)
     call groups(i)%ptr%scalepositions(oldbox, simbox)
  end do
  Vn = volume(simbox)
  
  !! Calculate potential energy in the scaled system.
  call total_energy(groups, simbox, pair_ias, single_ias, &
       totenew, err)
  
  if (err /= 0) then
     !! Scale particles back to old coordinates.
     do i = 1, size(groups)
        call groups(i)%ptr%scalepositions(simbox, oldbox)
     end do
     simbox = oldbox
  else
     boltzmannn = totenew + pressure * Vn - real(nparticles, dp) * &
          temperature * log(Vn)  
     boltzmanno = e_total + pressure * Vo - real(nparticles, dp) * &
          temperature * log(Vo)
     call acceptchange(boltzmanno, boltzmannn, temperature, genstate, &
          isaccepted)
     if (isaccepted) then
        e_total = totenew
     else 
        !! Scale particles back to old coordinates
        do i = 1, size(groups)
           call groups(i)%ptr%scalepositions(simbox, oldbox)
        end do
        simbox = oldbox
     end if
  end if

  if (present(n_accepted)) then
     if (isaccepted) then
        n_accepted = 1
     else
        n_accepted = 0
     end if
  end if

  if(present(n_trials)) n_trials = 1
end subroutine movevol

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
