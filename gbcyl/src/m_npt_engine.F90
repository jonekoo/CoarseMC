module m_npt_engine
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
  use m_particlegroup
  use m_nvt_engine, only: temperature, newmaxvalue
  implicit none  

  !> The simulation pressure. Only meaningful in a constant-pressure
  !! simulation.
  real(dp), save :: pressure = -1._dp
  
!> The maximum absolute change of volume in a trial volume update.
  real(dp), save :: maxscaling = 100._dp
  
  !> The types of trial volume updates.
  character(len = 200), dimension(:), allocatable, save :: scalingtypes
  
  !> True if the module is correctly initialized.
  logical :: is_initialized = .false.

  !> Counter for trial volume updates.
  integer, save :: nscalingtrials = 0
  
  !> The number of accepted volume moves
  integer, save :: nacceptedscalings = 0
  
  !> The desired acceptance ratio for trial volume updates.
  real(dp), save :: scalingratio
    
contains
  
  !> Initializes the module by getting the module parameters from the
  !! @p reader object.
  !!
  !! @param reader the object which gets parameters by name e.g. from an
  !!        input file.
  !! @param the_simbox the simulation box defining the boundary
  !!        conditions and volume.
  !! @param the_particles the particles to be updated with the MC moves.
  !!
  subroutine npt_engine_init(reader)
    type(parameterizer), intent(in) :: reader
    character(len = 200), save :: scalingtype = "z"
    call getparameter(reader, 'scaling_type', scalingtype)
    call parsescalingtype(scalingtype)
    call getparameter(reader, 'max_scaling', maxscaling)
    is_initialized = .true.
  end subroutine npt_engine_init
  
  subroutine npt_engine_init_json(json_val)
    type(json_value), pointer, intent(in) :: json_val
    allocate(scalingtypes(0))
    call get_parameter(json_val, 'scaling_types', scalingtypes) 
    call get_parameter(json_val, 'max_scaling', maxscaling, error_lb=0._dp)
    call get_parameter(json_val, 'pressure', pressure, error_lb=0._dp)
    call get_parameter(json_val, 'nscalingtrials', nscalingtrials, error_lb=0)
    call get_parameter(json_val, 'nacceptedscalings', nacceptedscalings, &
         error_ub=nscalingtrials, error_lb=0)
    call get_parameter(json_val, 'scaling_ratio', scalingratio, &
         error_lb=0._dp, error_ub=1._dp)
    is_initialized = .true.
  end subroutine npt_engine_init_json
  
  !> Finalizes the module state.
  subroutine npt_engine_finalize
    if (allocated(scalingtypes)) deallocate(scalingtypes)
    is_initialized = .false.
  end subroutine npt_engine_finalize
  
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
  subroutine npt_engine_to_json(json_val)
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
    call json_add(json_val, 'pressure', pressure)
    call json_add(json_val, 'scaling_ratio', scalingratio)
    call json_add(json_val, 'nscalingtrials', nscalingtrials)
    call json_add(json_val, 'nacceptedscalings', nacceptedscalings)
    if (nscalingtrials > 0) then
       call json_add(json_val, 'current_scaling_ratio', &
            real(nacceptedscalings, dp)/real(nscalingtrials, dp))
    else
       call json_add(json_val, 'current_scaling_ratio', 'nan')
    end if
  end subroutine npt_engine_to_json

  subroutine npt_engine_reset_counters()
      nscalingtrials = 0 
      nacceptedscalings = 0
  end subroutine npt_engine_reset_counters

  subroutine npt_update(simbox, groups, genstate, pair_ias, single_ias, &
       e_total)
    type(poly_box), intent(inout) :: simbox
    type(particlegroup_ptr), intent(inout) :: groups(:)
    type(rngstate), intent(inout) :: genstate
    type(pair_interaction_ptr), intent(in) :: pair_ias(:, :)
    type(single_interaction_ptr) :: single_ias(:)
    real(dp), intent(out) :: e_total
    integer :: n_trials, n_accepted
    call update_volume(simbox, groups, genstate, pair_ias,&
         single_ias, temperature, pressure, e_total, n_trials, &
         n_accepted)
    nscalingtrials = nscalingtrials + n_trials
    nacceptedscalings = nacceptedscalings + n_accepted
  end subroutine npt_update

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

  subroutine npt_engine_update_max_scaling()
    if (nscalingtrials > 0) then
       call set_maxscaling(newmaxvalue(nscalingtrials, nacceptedscalings, &
            scalingratio, get_maxscaling()))
    end if
  end subroutine npt_engine_update_max_scaling

end module m_npt_engine
