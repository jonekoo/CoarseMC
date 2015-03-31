module m_particle
  use nrtype
  use gen_lists
  use json_module
  use iso_fortran_env
  
  type particle
     type(list) :: interactions
     real(dp) :: position(3) = [0., 0., 0.]
   contains
     procedure :: add_interaction => add_gen_ia
     !   final :: particle_finalize
  end type particle

  interface particle
     module procedure from_json
  end interface particle
  
contains

  function from_json(json_val)
    type(particle) :: from_json
    type(json_value), intent(in), pointer :: json_val
    real(REAL64), allocatable :: temp(:)
    logical :: found
    call json_get(json_val, 'position', temp, found)
    from_json%position(1:3) = temp(1:3)
  end function from_json
  
  subroutine add_gen_ia(this, ia)
    class(particle), intent(inout) :: this
    class(*), intent(in) :: ia
    call this%interactions%add(ia)
  end subroutine add_gen_ia
  
  subroutine particle_finalize(this)
    type(particle), intent(inout) :: this
    call this%interactions%finalize()
  end subroutine particle_finalize
  
end module m_particle



module m_rod
  use m_particle
  use iso_fortran_env
  
  type, extends(particle) :: rod
     real(dp) :: orientation(3) = [0., 0., 1.]
   contains
     procedure :: add_interaction => add_rod_ia
  end type rod

  interface rod
     module procedure rod_from_json
  end interface rod
  
  
  type, abstract :: rod_interaction
   contains
     procedure(interact), deferred :: value
  end type rod_interaction
  
  interface
     subroutine interact(this, position, orientation, res, err)
       import dp
       import rod_interaction
       class(rod_interaction), intent(in) :: this
       real(dp), intent(in) :: position(3), orientation(3)
       real(dp), intent(out) :: res
       integer, intent(out) :: err
     end subroutine interact
  end interface

contains

  function rod_from_json(json_val)
    type(rod) :: rod_from_json
    type(json_value), pointer, intent(in) :: json_val
    real(REAL64), allocatable :: temp(:)
    logical :: found
    rod_from_json%particle = particle(json_val)
    call json_get(json_val, 'orientation', temp, found)
    rod_from_json%orientation(1:3) = temp(1:3)
  end function
  
  subroutine add_rod_ia(this, ia)
    class(rod), intent(inout) :: this
    class(*), intent(in) :: ia
    select type(ia)
    class is (rod_interaction)
       call this%particle%add_interaction(ia)
    class default
       write(error_unit, *) 'Error: Tried to add a non-compatible ' // &
            'interaction to a rod. Stopping.'
       stop
    end select
  end subroutine add_rod_ia
  
  subroutine energy(this, res, err)
    class(rod), intent(inout) :: this
    real(dp), intent(out) :: res
    integer, intent(out) :: err
    class(*), pointer :: ia
    real(dp) :: cur
    res = 0.
    call this%interactions%iter_restart()
    do while(.true.)
       call this%interactions%iter_next(ia)
       if (.not. associated(ia)) exit
       select type(ia)
       class is (rod_interaction)
          call ia%value(this%position, this%orientation, cur, err)
          if (err /= 0) exit
          res = res + cur
       end select
    end do
  end subroutine energy
  
end module m_rod



module m_particlegroup
  use nrtype
  use gen_lists
  use m_particle
  use json_module
  
  type, extends(particle) :: particlegroup
     character(kind=CK, len=:), allocatable :: name
     real(dp), allocatable :: positions(:, :)
   contains
     procedure :: to_json => particlegroup_to_json
  end type particlegroup
  
  interface particlegroup
     module procedure empty_group, particlegroup_from_json
  end interface particlegroup
  
contains 
  
  function empty_group(name)
    type(particlegroup) :: empty_group
    character(len=:), allocatable, intent(in) :: name
    empty_group%name = name(1:80)
  end function empty_group

  function particlegroup_from_json(json_val) result(pg)
    type(particlegroup) :: pg
    type(json_value), pointer, intent(in) :: json_val
    type(json_value), pointer :: child, gchild
    logical :: found
    real(REAL64), allocatable :: temp(:)
    character(kind=CK, len=:), allocatable :: temp_name
    call json_get(json_val, 'name', temp_name, found)
    if (.not. found) then
       write(error_unit, *) 'Warning: Creating a particlegroup without a name.'
    end if
    pg%name = temp_name
    call json_get(json_val, 'positions', child, found)
    if (.not. found) then
       write(error_unit, *) 'Warning: Creating an empty particlegroup.'
    else
       call get_3d_points_json(child, pg%positions)
    end if
  end function particlegroup_from_json

  subroutine get_3d_points_json(json_val, points)
    type(json_value), pointer, intent(in) :: json_val
    real(REAL64), allocatable, intent(out) :: points(:, :)
    type(json_value), pointer :: point_json
    integer :: i
    real(REAL64), allocatable :: temp(:)
    allocate(points(3, json_count(json_val)))
    do i = 1, json_count(json_val)
       call json_get_child(json_val, i, point_json)
       call json_get(point_json, temp)
       points(:, i) = temp(1:3)
    end do
  end subroutine get_3d_points_json

  subroutine particlegroup_to_json(this, json_val)
    class(particlegroup), intent(in) :: this
    type(json_value), pointer :: json_val
    character(kind=CK, len=:), allocatable :: temp_name
    logical :: found
    call json_create_object(json_val, '')
    call json_add(json_val, 'name', this%name)
    call add_3d_points_json(json_val, 'positions', this%positions)
  end subroutine particlegroup_to_json

  subroutine add_3d_points_json(json_val, name, points)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    real(REAL64), intent(in) :: points(:, :)
    type(json_value), pointer :: temp_array => null()
    call json_create_array(temp_array, name)
    do i = 1, size(points(1, :)) 
       call json_add(temp_array, '', points(:, i))
    end do
    call json_add(json_val, temp_array)
  end subroutine add_3d_points_json
  
end module m_particlegroup



module m_rodgroup
  use m_rod
  use m_particlegroup
  
  type, extends(particlegroup) :: rodgroup
     real(dp), allocatable :: orientations(:, :)
  end type rodgroup

  interface rodgroup
     module procedure rodgroup_json
  end interface rodgroup

contains

  function rodgroup_json(json_val)
    type(rodgroup) :: rodgroup_json
    type(json_value), pointer, intent(in) :: json_val
    type(json_value), pointer :: child
    logical :: found
    rodgroup_json%particlegroup = particlegroup(json_val)
    call json_get(json_val, 'orientations', child, found)
    call get_3d_points_json(child, rodgroup_json%orientations)
  end function rodgroup_json

  
end module m_rodgroup


module m_builder
  
  use m_particlegroup
  use iso_fortran_env
  use json_module
  use m_rod
  use gen_lists
  use m_rodgroup
  implicit none
  
contains 
  
  
  subroutine build(filename, groups)
    character(len=*), intent(in) :: filename
    !class(particle), allocatable, intent(out) :: groups(:)
    type(list) :: groups
    type(json_file) :: json
    logical :: found
    !integer :: i,j,k
    type(json_value), pointer :: p, child, pos, temp
    integer(INT32) :: i,j
    character(kind=CK, len=1) :: ic
    character(kind=CK, len=:), allocatable :: name
    character(kind=CK, len=:), allocatable :: type_str
    integer(INT32) :: group_size
    real(REAL64), allocatable :: positions(:, :), temp_pos(:)
    real(REAL64) :: x
    !real(REAL64), allocatable :: temp(:)
    class(particle), allocatable :: temp_particle
    call json_initialize()

    ! read particlegroups from file        
    call json%load_file(filename = filename)!'groups.json')
    call json%get('groups', p, found)
    !allocate(groups(json_count(p)))
    do i = 1, json_count(p)
       write(output_unit, *) 'Get particle from json:'
       call json_get_child(p, i, child)
       write(output_unit, *) 'Read particle from json:'
       call particle_factory(child, temp_particle)
       write(output_unit, *) 'Add particle to groups:'
       select type (temp_particle)
       class is (particlegroup)
          write(*, *) 'name, positions:', temp_particle%name, &
               temp_particle%positions
       class is (rodgroup)
          write(*, *) 'positions, orientations:', &
               temp_particle%positions, temp_particle%orientations
       end select
       call groups%add(temp_particle)
    end do
    
    ! clean up
    call json%destroy()
    if (json_failed()) stop 1
    
    ! build particlegrops
  end subroutine build


  subroutine particle_factory(json_val, p)
    type(json_value), pointer, intent(in) :: json_val
    class(particle), allocatable, intent(out) :: p
    character(kind=CK, len=:), allocatable :: type_str
    logical :: found
    call json_get(json_val, 'type', type_str, found)
    if (.not. found) then
       write(error_unit, *) 'Error: Could not read particle type!'
       stop 1
    end if
    type_str = adjustl(type_str)
    if (trim(type_str) == "particle") then
       allocate(p, source=particle(json_val))
    else if (trim(type_str) == "rod") then
       allocate(p, source=rod(json_val))
    else if (trim(type_str) == "particlegroup") then
       allocate(p, source=particlegroup(json_val))
    else if (trim(type_str) == "rodgroup") then
       allocate(p, source=rodgroup(json_val))
    end if
  end subroutine particle_factory

  subroutine serialize(groups, filename)
    type(list), intent(inout) :: groups
    character(len=*), intent(in) :: filename
    class(*), pointer :: res
    integer :: n = 0
    type(json_value), pointer :: json_val
    character(kind=CK, len=:), allocatable :: str
    call groups%iter_restart()
    do while(.true.)
       call groups%iter_next(res)
       if (.not. associated(res)) then
          write(*, *) 'end of list reached after', n, 'items'
          exit
       end if
       n = n + 1
       select type (res)
       class is (particlegroup)
          write(*, *) 'iterated ', res%name
          call res%to_json(json_val)
          !call json_print(json_val, filename)
          call json_print_to_string(json_val, str)
          write(*, *) str
       end select
    end do
  end subroutine serialize
  
end module m_builder



program test_composition
  use m_builder
  use gen_lists
  use json_module
  ! read interactions from file
  ! attach interactions to particlegroups
  !class(particle), allocatable :: groups(:)
  type(list) :: groups
  type(json_value), pointer :: group_json
  call build('groups.json', groups)
  call serialize(groups, 'output.json')
  
  ! simulate
end program test_composition

