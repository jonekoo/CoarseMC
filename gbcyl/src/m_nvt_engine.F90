!> Implements Metropolis Monte Carlo updates of the coordinates of the
!! particles and the simulation box dimensions. The module also
!! controls the parallel tempering updates. If compiled and run
!! with OpenMP, the domain decomposition algorithm is used for the
!! particle moves and energy calculations. 
module m_nvt_engine
  use num_kind, only: dp
  use iso_fortran_env, only: output_unit, error_unit
  use class_poly_box, only: poly_box
  use m_particle, only: particle, &
       pair_interaction, pair_interaction_ptr, &
       single_interaction, single_interaction_ptr, &
       particlearray_wrapper
  use utils, only: splitstr, join, acceptchange
  !$ use omp_lib
  use class_simplelist, only: simplelist, new_simplelist, simplelist_update, &
       simplelist_nbr_cells, flat_index, simplelist_deallocate, &
       simplelist_nbrmask, simplelist_cell_nbrmask
  include 'rng.inc'
  use json_module, only: json_get, json_value, json_add
  use m_json_wrapper, only: get_parameter
  use m_particlegroup, only: particlegroup, particlegroup_ptr
  use particle_mover, only: setmaxmoves, getmaxmoves
  implicit none  

  type, extends(particlearray_wrapper) :: domain
    integer :: n_cell = 0
  contains
    procedure :: domain_assign
    generic :: assignment(=) => domain_assign
    final :: domain_delete
  end type

  !> Counter for trial particle moves.
  integer, save :: nmovetrials = 0
  
  !> The number of accepted particle moves
  integer, save :: nacceptedmoves = 0
  
  !> The desired acceptance ratio for trial particle moves.
  real(dp), save :: moveratio
  
  !> The simulation temperature.
  real(dp), save :: temperature = -1._dp

contains

  subroutine nvt_engine_init_json(json_val)
    type(json_value), intent(in), pointer :: json_val
    call get_parameter(json_val, 'move_ratio', moveratio, error_lb=0._dp, &
         error_ub=1._dp)
    call get_parameter(json_val, 'temperature', temperature, error_lb=0._dp)
    call get_parameter(json_val, 'nmovetrials', nmovetrials, error_lb=0)
    call get_parameter(json_val, 'nacceptedmoves', nacceptedmoves, &
         error_ub=nmovetrials, error_lb=0)
  end subroutine nvt_engine_init_json

  subroutine nvt_engine_to_json(json_val)
    type(json_value), pointer, intent(inout) :: json_val
    call json_add(json_val, 'temperature', temperature)
    call json_add(json_val, 'move_ratio', moveratio)
    call json_add(json_val, 'nmovetrials', nmovetrials)
    call json_add(json_val, 'nacceptedmoves', nacceptedmoves)
    if (nmovetrials > 0) then
       call json_add(json_val, 'current_move_ratio', &
            real(nacceptedmoves, dp)/real(nmovetrials, dp))
    else
       call json_add(json_val, 'current_move_ratio', 'nan')
    end if
  end subroutine nvt_engine_to_json

  
  subroutine domain_assign(this, src)
    class(domain), intent(inout) :: this
    type(domain), intent(in) :: src
    this%particlearray_wrapper = src%particlearray_wrapper
    this%n_cell = src%n_cell
  end subroutine domain_assign

  subroutine domain_delete(this)
    type(domain), intent(inout) :: this
    this%n_cell = 0
  end subroutine domain_delete


  subroutine nvt_update(groups, genstates, simbox, pair_ias, single_ias, dE)
    type(particlegroup_ptr), intent(inout) :: groups(:)
    type(poly_box), intent(in) :: simbox
    type(rngstate), intent(inout) :: genstates(0:)
    type(pair_interaction_ptr), intent(in) :: pair_ias(:, :)
    type(single_interaction_ptr) :: single_ias(:)
    real(dp), intent(out) :: dE
    integer :: n_trials, n_accepted
    call make_particle_moves(groups, genstates, simbox, temperature, & 
         pair_ias, single_ias, dE, n_trials, n_accepted)
    nmovetrials = nmovetrials + n_trials
    nacceptedmoves = nacceptedmoves + n_accepted
  end subroutine nvt_update

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
#ifdef DEBUG
       if (.not. groups(i_group)%ptr%check_particles(simbox)) then
          write(error_unit, *) 'groups(', i_group, &
               '): some particles are not in the box before moves!'
          stop
       end if
#endif
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
                         ds(i_group) = create_domain(groups(i_group)%ptr, &
                              simbox, jx, jy, jz)
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
                         !! dummy argument that has intent(out) or 
                         !! intent(inout)."
                         do i = 1, ds(i_group)%n_cell
                            call groups(i_group)%ptr%particles(&
                                 groups(i_group)%ptr%sl%indices(i, jx, jy, jz)&
                                 )%downcast_assign(ds(i_group)%arr(i))
                         end do
                      end do
                      do i_group = 1, size(groups)
                         call domain_delete(ds(i_group))
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
    !$OMP END PARALLEL

    do i_group = 1, size(groups)
#ifdef DEBUG
       if (.not. groups(i_group)%ptr%check_particles(simbox)) then
          write(error_unit, *) 'groups(', i_group, &
               '): some particles are not in the box!'
          stop
       end if
#endif
       call simplelist_update(groups(i_group)%ptr%sl, simbox, &
            groups(i_group)%ptr%particles)
    end do
  end subroutine make_particle_moves
  
  
  subroutine nvt_engine_reset_counters()
    nacceptedmoves = 0
    nmovetrials = 0
  end subroutine nvt_engine_reset_counters
  
  
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
    !! array class(particle) so we'll use a workaround.
    !! Put particles in jx, jy, jz first in temp_particles:
    allocate(helper, source=[this%sl%indices(1:d%n_cell, jx, jy, jz), &
         pack([(i, i = 1, size(this%particles))], nbr_mask)])
    
    allocate(d%arr(n_local), mold=this%particles(1:n_local))
    do i = 1, n_local
       call d%arr(i)%downcast_assign(this%particles(helper(i)))
    end do
    
    allocate(d%mask(n_local), source=.true.)
  end function create_domain
  

  subroutine delete_domain(this)
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
    class(particle), allocatable :: newparticle
    integer :: err
    real(dp) :: enew
    real(dp) :: eold
    logical :: isaccepted
    n_accepted = 0
    n_trials = 0
    dE = 0.
    if (domains(i_d)%n_cell > 0) allocate(newparticle, &
         source=domains(i_d)%arr(1))
    do j = 1, domains(i_d)%n_cell
       domains(i_d)%mask(j) = .false.
#ifdef DEBUG
       call newparticle%downcast_assign(domains(i_d)%arr(j), err)
       if (err /= 0) then
          write(error_unit, *) 'newparticle%downcast_assign: err=', err
          stop 'domain_move unable to continue'
       end if
#else
       call newparticle%downcast_assign(domains(i_d)%arr(j))
#endif
       call newparticle%move(genstates(1))
       call newparticle%set_position(simbox%minimage(newparticle%x, &
            newparticle%y, newparticle%z))
       call newparticle%energy(domains, pair_ias(:, i_d), simbox, &
            single_ias(i_d)%ptr, enew, err)
       if(err == 0) then 
          call domains(i_d)%arr(j)%energy(domains, pair_ias(:, i_d), simbox, &
               single_ias(i_d)%ptr, eold, err)
          if (err /= 0) then
             write(error_unit, *) 'ERROR: err=', err, ' domain=', i_d, &
                  ' old particle ', j
             if (err == 1) write(error_unit, *) &
                  'Particle energy calculation resulted in overlap before move.'
             stop 'Stopped by domain_move.'
          end if
          call acceptchange(eold, enew, temperature, genstates(1), isaccepted)
          if(isaccepted) then
#ifdef DEBUG
             call domains(i_d)%arr(j)%downcast_assign(newparticle, err)
             if (err /= 0) then
                write(error_unit, *) 'domains(', i_d, &
                     '%arr(', j, ')%downcast_assign: err = ', err
                stop 'domain_move unable to continue'
             end if
#else
             call domains(i_d)%arr(j)%downcast_assign(newparticle)
#endif
             dE = dE + enew - eold
             n_accepted = n_accepted + 1
          end if
       end if
       domains(i_d)%mask(j) = .true.
    end do
    n_trials = n_trials + domains(i_d)%n_cell
  end subroutine domain_move

  subroutine nvt_engine_update_max_moves()
    real(dp) :: newdthetamax
    real(dp) :: newdximax
    if (nmovetrials > 0) then
       call getmaxmoves(newdximax, newdthetamax)
       !! Adjust translation
       newdximax = newmaxvalue(nmovetrials, nacceptedmoves, moveratio, &
            newdximax)
       !! Update the minimum cell side length of the cell list because the
       !! maximum translation has changed: 
     
       !! This should adjust rotations < pi/2 to move the particle end as
       !! much as a random translation. 4.4 is the assumed molecule 
       !! length-to-breadth ratio.
       newdthetamax = 2._dp * asin(newdximax/4.4_dp) 
       call setmaxmoves(newdximax, newdthetamax)
    end if
  end subroutine nvt_engine_update_max_moves
  
  !> Returns a new trial move parameter value calculated from the desired
  !! acceptance ratio.
  !! 
  !! @param ntrials the total number of trials of the kind of move in
  !!        question.
  !! @param naccepted number of accepted trials.
  !! @param desiredratio the desired acceptance ratio for trial moves.
  !! @param oldvalue the old value of the parameter setting the maximum
  !!        size for the trial move in question.
  !!
  function newmaxvalue(ntrials, naccepted, desiredratio, oldvalue) &
       result(newvalue)
    integer, intent(in) :: ntrials
    integer, intent(in) :: naccepted
    real(dp), intent(in) :: desiredratio
    real(dp), intent(in) :: oldvalue
    real(dp) :: newvalue
    real(dp), parameter :: multiplier = 1.05_dp
    if (real(naccepted, dp) / real(ntrials, dp) > desiredratio) then
       newvalue = oldvalue * multiplier
    else 
       newvalue = oldvalue / multiplier
    end if
  end function newmaxvalue
  
end module m_nvt_engine
