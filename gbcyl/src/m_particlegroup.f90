!> Implements Metropolis Monte Carlo updates of the coordinates of the
!! particles and the simulation box dimensions. The module also
!! controls the parallel tempering updates. If compiled and run
!! with OpenMP, the domain decomposition algorithm is used for the
!! particle moves and energy calculations. 
module m_particlegroup
  use nrtype
  use class_poly_box
  use particle
  use class_parameterizer
  use class_parameter_writer
  use genvoltrial
  use utils 
  !$ use omp_lib
  use class_simplelist
  include 'rng.inc'
  implicit none  

  !> The maximum absolute change of volume in a trial volume update.
  real(dp), save :: maxscaling = 100._dp

  !> The total energy.
  real(dp), save :: etotal = 0._dp

  !> The types of trial volume updates.
  character(len = 200), dimension(:), allocatable, save :: scalingtypes

  !> True if the module is correctly initialized.
  logical :: is_initialized = .false.

  type particlegroup
     type(particledat), allocatable :: particles(:)
     type(simplelist) :: sl
   contains
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
  end type

contains

  function create_particlegroup(simbox, particles, min_cell_length, &
       min_boundary_width) result(group)
    type(poly_box), intent(in) :: simbox
    type(particledat), intent(in) :: particles(:)
    real(dp), intent(in) :: min_cell_length, min_boundary_width
    type(particlegroup) :: group
    allocate(group%particles(size(particles)), source=particles)
    !$ if (.true.) then
    !$ call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
    !$& min_boundary_width, is_x_even = isxperiodic(simbox), &
    !$& is_y_even = isyperiodic(simbox), is_z_even = iszperiodic(simbox))
    !$ else 
    call new_simplelist(group%sl, simbox, group%particles, min_cell_length, &
         min_boundary_width)
    !$ end if
  end function create_particlegroup

  impure elemental subroutine particlegroup_finalize(group)
    type(particlegroup), intent(inout) :: group
    call simplelist_deallocate(group%sl)
    if(allocated(group%particles)) deallocate(group%particles)
  end subroutine particlegroup_finalize
  
  impure elemental subroutine domain_assign(dest, src)
    class(domain), intent(inout) :: dest
    type(domain), intent(in) :: src
    dest%particlearray_wrapper = src%particlearray_wrapper
    dest%n_cell = src%n_cell
  end subroutine domain_assign

  
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

!> Finalizes the module state.
subroutine mcsweep_finalize
  if (allocated(scalingtypes)) deallocate(scalingtypes)
  is_initialized = .false.
end subroutine mcsweep_finalize

!> Returns the total energy of the system.
real(dp) function get_total_energy()
  if (.not. is_initialized) then
     stop 'mc_sweep:get_total_energy: Error: module not initialized!'
  end if
  get_total_energy = etotal
end function get_total_energy

!> Checks and parses the @p scalingtype string to the scalingtypes
!! array. 
subroutine parsescalingtype(scalingtype)
  character(len=*), intent(in) :: scalingtype
  logical :: isok
  isok = (0 == verify(trim(adjustl(scalingtype)), 'xyz,') .and. &
       0 /= len_trim(adjustl(scalingtype)))
  if (.not. isok) then
     write(*, *) 'Error parsing parameter scaling_type, stopping!'
     stop
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
  call writeparameter(writer, 'total_energy', etotal)
end subroutine mc_sweep_writeparameters

!> Schedules parallel moves of particles using OpenMP with a domain 
!! decomposition algorithm. 
!!
!! @see e.g. G. Heffelfinger and M. Lewitt. J. Comp. Chem., 17(2):250–265,
!! 1996.
!! 
!! @param simbox is the simulation box in which the particles reside.
!! @param particles is the array of particles to move.
!! @param genstates random number generator states. Each thread needs one.
!! @param sl the cell list presenting the decomposition.
!! 
subroutine make_particle_moves(groups, genstates, simbox, temperature, &
     pair_ia, subr_single_energy, dE, n_trials, n_accepted)
  implicit none
  class(particlegroup_ptr), intent(inout) :: groups(:)
  type(poly_box), intent(in) :: simbox
  type(rngstate), intent(inout) :: genstates(0:)
  class(pair_interaction), intent(in) :: pair_ia
  procedure(single_energy) :: subr_single_energy
  real(dp), intent(in) :: temperature
  real(dp), intent(out) :: dE
  integer, intent(out) :: n_accepted, n_trials
  
  !$ integer :: n_threads
  integer :: thread_id
  !! You may be tempted to make the allocatable arrays automatic, but there's
  !! no performance gained and depending on the system large automatic arrays
  !! may cause a stack overflow.
  real(dp) :: dE_d
  integer :: ix, iy, iz, jx, jy, jz, n_trials_d, n_accepted_d, i_group
  type(domain), allocatable :: ds(:)
  do i_group = 1, size(groups)
     call update(groups(i_group)%ptr%sl, simbox, groups(i_group)%ptr%particles)
  end do
  thread_id = 0
  !$ n_threads = 1
  dE = 0._dp
  dE_d = 0._dp
  n_accepted = 0
  n_trials = 0
  !! Loop over cells. This can be thought of as looping through a 
  !! 2 x 2 x 2 cube of cells.
  !$OMP PARALLEL shared(groups, simbox, genstates, pair_ia)& 
  !$OMP& private(thread_id, n_threads, n_trials_d, n_accepted_d, ds)&
  !$OMP& reduction(+:dE, n_accepted, n_trials) 
  !$ thread_id = omp_get_thread_num()
  allocate(ds(size(groups)))
  do iz=0, min(1, groups(1)%ptr%sl%nz-1)
     do iy=0, min(1, groups(1)%ptr%sl%ny-1)
        do ix=0, min(1, groups(1)%ptr%sl%nx-1)
           !$OMP DO collapse(3) private(i_group, dE_d)&
           !$OMP& schedule(dynamic)
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
                       call domain_move(ds(i_group), &
                            genstates(thread_id:thread_id), simbox, &
                            temperature, &
                            [ds(:i_group - 1), ds(i_group + 1:)], pair_ia, &
                            subr_single_energy, dE_d, n_trials=n_trials_d, &
                            n_accepted=n_accepted_d)
                       dE = dE + dE_d
                       n_trials = n_trials + n_trials_d
                       n_accepted = n_accepted + n_accepted_d
                    end do
                    !! Synchronize
                    !! :TODO: Generalize to a sync domain operation.
                    do i_group = 1, size(groups)   
                       groups(i_group)%ptr%particles(&
                            groups(i_group)%ptr%sl%indices(&
                            1:ds(i_group)%n_cell, jx, jy, jz)) = &
                            ds(i_group)%arr(1:ds(i_group)%n_cell)
                    end do
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
  !if (allocated(temp_particles)) deallocate(temp_particles)
  !if (allocated(nbr_mask)) deallocate(nbr_mask)
  !$OMP END PARALLEL
  etotal = etotal + dE
  do i_group = 1, size(groups)
     call update(groups(i_group)%ptr%sl, simbox, groups(i_group)%ptr%particles)
  end do
!  call update(this%sl, simbox, this%particles)
end subroutine make_particle_moves


function create_domain(this, simbox, jx, jy, jz) result(d)
  type(particlegroup), intent(in) :: this
  type(poly_box), intent(in) :: simbox
  integer, intent(in) :: jx, jy, jz
  type(domain) :: d
  logical, allocatable :: nbr_mask(:)
  integer :: n_local
  allocate(nbr_mask(size(this%particles)), source=.false.)
  d%n_cell = this%sl%counts(jx, jy, jz)
  call simplelist_nbrmask(this%sl, simbox, jx, jy, jz, &
       nbr_mask)
  n_local = count(nbr_mask)
  
  !! Remove particles in jx, jy, jz from mask:
  nbr_mask(this%sl%indices(1:d%n_cell, jx, jy, jz)) = .false.
  !! Put particles in jx, jy, jz first in temp_particles:
  allocate(d%arr(n_local), source=[this%particles(&
       this%sl%indices(1:d%n_cell, jx, jy, jz)), &
       pack(this%particles, nbr_mask)])
  allocate(d%mask(size(d%arr)), source=.true.)
end function create_domain

impure elemental subroutine delete_domain(this)
  class(domain), intent(inout) :: this
  this%n_cell = 0
end subroutine delete_domain

subroutine domain_move(this, genstates, simbox, temperature, nbrs, &
     pair_ia, subr_single_energy, dE, n_trials, n_accepted)
  type(domain), intent(inout) :: this
  type(rngstate), intent(inout) :: genstates(:)
  type(poly_box), intent(in) :: simbox
  real(dp), intent(in) :: temperature
  type(domain), intent(in) :: nbrs(:)
  class(pair_interaction), intent(in) :: pair_ia
  procedure(single_energy) :: subr_single_energy
  real(dp), intent(out) :: dE
  integer, intent(out), optional :: n_trials, n_accepted
  integer :: j
  real(dp) :: dE_j
  class(particledat), allocatable :: newparticle
  class(particledat), allocatable :: oldparticle

  integer :: err
  real(dp) :: enew
  real(dp) :: eold
  logical :: isaccepted
  class(particlearray_wrapper), allocatable :: temp(:)
  n_accepted = 0
  n_trials = 0
  !! GFortran 6.0.0 (BETA) needs the allocation for the code to work.
  !! ifort 16.0.0 works without
  allocate(temp(size(nbrs) + 1), source=[this, nbrs])
  do j = 1, this%n_cell
     temp(1)%mask(j) = .false.
     enew = 0._dp
     eold = 0._dp
     dE_j = 0._dp
     isaccepted = .false.
     if (allocated(newparticle)) deallocate(newparticle)
     allocate(newparticle, source=this%arr(j))
     call move(newparticle, genstates(1))
     call setposition(newparticle, minimage(simbox, position(newparticle)))
     if (allocated(oldparticle)) deallocate(oldparticle)
     allocate(oldparticle, source=this%arr(j))
     this%arr(j) = newparticle
 
     call newparticle%energy(temp, pair_ia, simbox, &
          subr_single_energy, enew, err)
   
     this%arr(j) = oldparticle
     if(err == 0) then 
        call oldparticle%energy(temp, pair_ia, simbox, &
             subr_single_energy, enew, err)
         call acceptchange(eold, enew, temperature, genstates(1), isaccepted)
        if(isaccepted) then
           this%arr(j) = newparticle
           dE_j = enew - eold
        end if
     end if
     
     if (isaccepted) n_accepted = n_accepted + 1
     n_trials = n_trials + 1
     
     dE = dE + dE_j
     temp(1)%mask(j) = .true.
  end do
end subroutine domain_move

subroutine update_volume(simbox, groups, genstate, pair_ia, subr_single_energy,&
     temperature, pressure, n_trials, n_accepted)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup_ptr), intent(inout) :: groups(:)
  type(rngstate), intent(inout) :: genstate
  class(pair_interaction), intent(in) :: pair_ia
  procedure(single_energy) :: subr_single_energy
  real(dp), intent(in) :: temperature, pressure
  integer, intent(out), optional :: n_trials, n_accepted
  integer :: i
  do i = 1, size(scalingtypes)
     call movevol(simbox, groups, scalingtypes(i), genstate, pair_ia, &
          subr_single_energy, temperature, pressure, n_trials, n_accepted)
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
subroutine movevol(simbox, groups, scalingtype, genstate, pair_ia, &
     subr_single_energy, temperature, pressure, n_trials, n_accepted)
  type(poly_box), intent(inout) :: simbox
  type(particlegroup_ptr), intent(inout) :: groups(:)
  character(len=*), intent(in) :: scalingtype
  type(rngstate), intent(inout) :: genstate
  class(pair_interaction), intent(in) :: pair_ia
  procedure(single_energy) :: subr_single_energy
  real(dp), intent(in) :: temperature, pressure
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
  
  !! It seems that total energy may drift if it is not updated here:
  call total_energy(groups, simbox, pair_ia, subr_single_energy, &
       etotal, err)
  if (err /= 0) stop 'movevol: overlap in old configuration! '//&
       'Should never happen!'
  
  !! Scale coordinates and the simulations box
  scaling = genvoltrial_scale(simbox, maxscaling, genstate, &
       trim(adjustl(scalingtype)))
  call groups(1)%ptr%scalepositions(oldbox, simbox)
  Vn = volume(simbox)
  
  !! Calculate potential energy in the scaled system.
  call total_energy(groups, simbox, pair_ia, subr_single_energy, &
       totenew, err)
  
  if (err /= 0) then
     !! Scale particles back to old coordinates.
     call groups(1)%ptr%scalepositions(simbox, oldbox)
     simbox = oldbox
  else
     boltzmannn = totenew + pressure * Vn - real(nparticles, dp) * &
          temperature * log(Vn)  
     boltzmanno = etotal + pressure * Vo - real(nparticles, dp) * &
          temperature * log(Vo)
     call acceptchange(boltzmanno, boltzmannn, temperature, genstate, &
          isaccepted)
     if (isaccepted) then
        etotal = totenew
     else 
        !! Scale particles back to old coordinates
        call groups(1)%ptr%scalepositions(simbox, oldbox)
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
subroutine total_energy(groups, simbox, pair_ia, subr_single_energy, energy, &
     err)
  type(particlegroup_ptr), intent(inout) :: groups(:)
  type(poly_box), intent(in) :: simbox
  class(pair_interaction), intent(in) :: pair_ia
  procedure(single_energy) :: subr_single_energy
  real(dp), intent(out) :: energy
  integer, intent(out) :: err
  integer :: i, j, ix, iy, iz, i_group
  logical :: mask(size(groups(1)%ptr%particles))
  real(dp) :: energy_j
  integer :: err_j
  integer :: helper(size(groups(1)%ptr%particles))
  integer, allocatable :: temp_helper(:)
  integer :: n_mask
  integer :: temp_j
  type(particledat), allocatable :: temp_particles(:)
  integer :: nbr_cells(3, 27), n_nbr_cells
  integer :: j_group
  helper = (/(i, i=1, size(groups(1)%ptr%particles))/) !! ifort vectorizes
  energy = 0._dp
  err = 0
  
  !$OMP PARALLEL default(shared) reduction(+:energy) reduction(+:err)& 
  !$OMP& private(energy_j, err_j, i, j, mask, temp_particles, temp_j, temp_helper, n_mask, nbr_cells, n_nbr_cells)
  allocate(temp_particles(size(groups(1)%ptr%particles)), &
       temp_helper(size(groups(1)%ptr%particles))) 
  !$OMP DO collapse(3) schedule(dynamic)
  do ix = 0, groups(1)%ptr%sl%nx - 1 
     do iy = 0, groups(1)%ptr%sl%ny - 1
        do iz = 0, groups(1)%ptr%sl%nz - 1

           do i_group = 1, size(groups)
              !! 1. compute inside ix, iy, iz in i_group
              call cell_energy(groups(i_group)%ptr, ix, iy, iz, simbox, &
                   pair_ia, subr_single_energy, energy_j, err)
              if (err /= 0) exit
              energy = energy + energy_j
           end do
           
           do i_group = 1, size(groups) - 1
              !! 2. compute with ix, ix, y in j_group > i_group
              do j_group = i_group + 1, size(groups)
                 call cell_pair_energy(groups(i_group)%ptr, ix, iy, iz, &
                      groups(j_group)%ptr, ix, iy, iz, simbox, pair_ia, &
                      energy_j, err)
                 if (err /= 0) exit
                 energy = energy + energy_j
              end do
              if (err /= 0) exit
           end do
           
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
                       call cell_pair_energy(groups(i_group)%ptr, ix, iy, iz, &
                            groups(j_group)%ptr, nbr_cells(1, i), &
                            nbr_cells(2, i), nbr_cells(3, i), simbox, pair_ia,&
                            energy_j, err)
                       if (err /= 0) exit
                       energy = energy + energy_j
                    end do
                 end if
                 if (err /= 0) exit
              end do
              if (err /= 0) exit
           end do
        end do
     end do
  end do
  !$OMP END DO
  if (allocated(temp_particles)) deallocate(temp_particles)
  if (allocated(temp_helper)) deallocate(temp_helper)
  !$OMP END PARALLEL  
end subroutine total_energy

subroutine cell_energy(this, ix, iy, iz, simbox, pair_ia, subr_single_energy, &
     energy, err)
  type(particlegroup), intent(in) :: this
  integer, intent(in) :: ix, iy, iz
  type(poly_box), intent(in) :: simbox
  class(pair_interaction), intent(in) :: pair_ia
  procedure(single_energy) :: subr_single_energy
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
               [particlej%x - particlei%x, particlej%y - particlei%y,&
               particlej%z - particlei%z])
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

  do i = 1, this%sl%counts(ix, iy, iz)
     associate(particlei => this%particles(this%sl%indices(i, ix, iy, iz)))
       call subr_single_energy(particlei, simbox, energy_ij, err)
       if (err /= 0) exit
       energy = energy + energy_ij
     end associate
  end do
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
            rij = minimage(simbox, [particlej%x - particlei%x, &
                 particlej%y - particlei%y, particlej%z - particlei%z])
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
  call update(this%sl, newbox, this%particles)
end subroutine scalepositions

end module m_particlegroup
