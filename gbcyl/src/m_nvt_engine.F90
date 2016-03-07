!> Implements the domain decomposition algorithm for moving particles.
!! The algorithm relies on the cell list implementation in
!! class_simplelist.
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
  use mt_stream, only: rng=>genrand_double1_s, rngstate=>mt_state
  use json_module, only: json_get, json_value, json_add
  use m_json_wrapper, only: get_parameter
  use m_particlegroup, only: particlegroup, particlegroup_ptr
  use particle_mover, only: setmaxmoves, getmaxmoves
  implicit none  

  !> Stores particles in one cell and its neighbour cells.
  type, extends(particlearray_wrapper) :: domain
     integer :: n_cell = 0
       !! Number of particles in the cell.
  contains
    procedure :: domain_assign
    generic :: assignment(=) => domain_assign
    procedure :: delete => domain_manual_delete
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

  !> Initializes the module by reading variables from JSON. Must be
  !! called before any other routine in this module.
  !!
  !! @param json_val contains the JSON.
  !! 
  subroutine nvt_engine_init_json(json_val)
    type(json_value), intent(in), pointer :: json_val
    call get_parameter(json_val, 'move_ratio', moveratio, error_lb=0._dp, &
         error_ub=1._dp)
    call get_parameter(json_val, 'temperature', temperature, error_lb=0._dp)
    call get_parameter(json_val, 'nmovetrials', nmovetrials, error_lb=0)
    call get_parameter(json_val, 'nacceptedmoves', nacceptedmoves, &
         error_ub=nmovetrials, error_lb=0)
  end subroutine nvt_engine_init_json

  !> Serializes the module state to JSON. The JSON is added to
  !! @p json_val.
  !! 
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

  !> Assignment operator implementation.
  subroutine domain_assign(this, src)
    class(domain), intent(inout) :: this
    type(domain), intent(in) :: src
    this%particlearray_wrapper = src%particlearray_wrapper
    this%n_cell = src%n_cell
  end subroutine domain_assign

  !> Final procedure.
  subroutine domain_delete(this)
    type(domain), intent(inout) :: this
    this%n_cell = 0
  end subroutine domain_delete

  !> Manual destructor when the final procedure is not applicable.
  subroutine domain_manual_delete(this)
    class(domain), intent(inout) :: this
    call this%particlearray_wrapper%delete()
    this%n_cell = 0
  end subroutine domain_manual_delete

  !> A wrapper around make_particle_moves.
  !!
  !! @param groups the particle groups to be updated.
  !! @param genstates the random number generator states for each
  !!        thread.
  !! @param simbox the simulation box.
  !! @param pair_ias the matrix of pair interactions.
  !! @param single_ias the array of single interactions.
  !! @param dE the change in energy after all trial moves have been
  !!        completed.
  !!
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
  !! @param groups the particle groups to be updated.
  !! @param genstates the random number generator states for each
  !!        thread.
  !! @param simbox the simulation box.
  !! @param pair_ias the matrix of pair interactions.
  !! @param single_ias the array of single interactions.
  !! @param dE the change in energy after all trial moves have been
  !!        completed.
  !! @param n_trials the number of trial moves attempted.
  !! @param n_accepted the number of accepted moves.
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
    integer :: ix, iy, iz, jx, jy, jz, n_trials_d, n_accepted_d, i_group, i,&
    	    n_max
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
    !$OMP& i_group, n_max)&
    !$OMP& reduction(+:dE, n_accepted, n_trials) 
    !$ thread_id = omp_get_thread_num()
    allocate(ds(size(groups)))
    do i_group = 1, size(groups)
      n_max = min(maxval(groups(i_group)%ptr%sl%counts) * 27, &
        size(groups(i_group)%ptr%particles))
      allocate(ds(i_group)%arr(n_max), &
        source=groups(i_group)%ptr%particles)
      allocate(ds(i_group)%mask(n_max), source = .false.)
      ds(i_group)%n_cell = 0
      ds(i_group)%n = 0
    end do  	 
    do iz=0, min(1, groups(1)%ptr%sl%nz-1)
       do iy=0, min(1, groups(1)%ptr%sl%ny-1)
          do ix=0, min(1, groups(1)%ptr%sl%nx-1)
             !$OMP DO collapse(3) schedule(dynamic)
             do jz = iz, groups(1)%ptr%sl%nz - 1, 2
                do jy = iy, groups(1)%ptr%sl%ny - 1, 2
                   do jx = ix, groups(1)%ptr%sl%nx - 1, 2
                      
                      !! Collect temp_particles for all groups
                      do i_group = 1, size(groups)
			 call set_domain(groups(i_group)%ptr, &
			       simbox, jx, jy, jz, ds(i_group))
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
                         do i = 1, ds(i_group)%n_cell
                            call groups(i_group)%ptr%particles(&
                                 groups(i_group)%ptr%sl%indices(i, jx, jy, jz)&
                                 )%downcast_assign(ds(i_group)%arr(i))
                         end do
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
    do i_group = 1, size(groups)
       call ds(i_group)%delete()
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
  
  !> Resets the counters used when adjusting moves.
  subroutine nvt_engine_reset_counters()
    nacceptedmoves = 0
    nmovetrials = 0
  end subroutine nvt_engine_reset_counters

  !> Assigns the cell @p jx @p jy @p jz to the domain @p d.
  !!
  !! @param this the particlegroup of the domain.
  !! @param simbox the simulation box.
  !! @param jx, jy, jz the cell index.
  !! @param d the domain with particles from cell jx, jy, jz at return.
  !!
  subroutine set_domain(this, simbox, jx, jy, jz, d)
    type(particlegroup), intent(in) :: this
    type(poly_box), intent(in) :: simbox
    integer, intent(in) :: jx, jy, jz
    type(domain), intent(inout) :: d
    logical, allocatable :: nbr_mask(:)
    integer :: i
    integer, allocatable :: helper(:)
    d%mask = .false.
    allocate(nbr_mask(size(this%particles)), source=.false.)
    d%n_cell = this%sl%counts(jx, jy, jz)
    call simplelist_cell_nbrmask(this%sl, simbox, jx, jy, jz, nbr_mask)
    d%n = count(nbr_mask)
    !! Remove particles in jx, jy, jz from mask:
    nbr_mask(this%sl%indices(1:d%n_cell, jx, jy, jz)) = .false.
    !! GFortran 6.0.0 is not fine with pack from a polymorphic
    !! array class(particle) so we'll use a workaround.
    allocate(helper, source=[this%sl%indices(1:d%n_cell, jx, jy, jz), &
         pack([(i, i = 1, size(this%particles))], nbr_mask)])
    if (allocated(nbr_mask)) deallocate(nbr_mask)
    do i = 1, d%n
       call d%arr(i)%downcast_assign(this%particles(helper(i)))
    end do
    d%mask(1:d%n) = .true.
  end subroutine set_domain  

  
  !> Moves particle in the cell of the domain @p domains(i_d). Used in
  !! the domain decomposition algorithm.
  !!
  !! @param domains are the subsystems of particles corresponding to
  !! 	    groups given to the make_particle_moves algorithm.
  !! @param i_d is the index of the domain to be moved in @p domains.
  !! @param genstates the random number generator states.
  !! @param simbox the simulation box.
  !! @param temperature the simulation temperature.
  !! @param pair_ias the pair interactions concerning domains(i_d)
  !! @param single_ias the single-particle interaction concerning 
  !! 	    domains(i_d). 
  !! @param dE the total change in energy due to the moves.
  !! @param n_trials the number of trial moves that were performed.
  !! @param n_accepted the number of accepted moves of n_trials.
  !!
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

  !> Adjusts the maximum moves of the particles.
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
