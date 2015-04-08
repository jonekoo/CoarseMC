module m_particle
  use nrtype
  use json_module
  use iso_fortran_env
  
  type particle
     real(dp) :: position(3) = [0., 0., 0.]
     contains
       procedure :: energy => particle_energy
  end type particle

  type, abstract :: particle_interaction
   contains
     procedure(interact), deferred :: value
  end type particle_interaction

  interface particle
     module procedure from_json
  end interface particle

  interface 
     pure subroutine interact(this, position, res, err)
       import particle_interaction
       import REAL64
       class(particle_interaction), intent(inout) :: this
       real(REAL64), intent(in) :: position(3)
       real(REAL64), intent(out) :: res
       integer, intent(out) :: err
     end subroutine interact
  end interface
  
contains

  function from_json(json_val)
    type(particle) :: from_json
    type(json_value), intent(in), pointer :: json_val
    real(REAL64), allocatable :: temp(:)
    logical :: found
    call json_get(json_val, 'position', temp, found)
    from_json%position(1:3) = temp(1:3)
  end function from_json
  
  elemental subroutine particle_energy(this, ia, res, err)
    class(particle), intent(in) :: this
    class(*), intent(inout) :: ia
    real(REAL64), intent(out) :: res
    integer, intent(out) :: err
    select type (ia)
    class is (particle_interaction)
      call ia%value(this%position, res, err)
    end select
  end subroutine particle_energy

end module m_particle



module m_rod
  use m_particle
  use iso_fortran_env

  type, extends(particle) :: rod
     real(dp) :: orientation(3) = [0., 0., 1.]
  end type rod

  interface rod
     module procedure rod_from_json
  end interface rod
    
  type, abstract :: rod_interaction
   contains
     procedure(interact_with_rod), deferred :: value
  end type rod_interaction
  
  interface
     pure subroutine interact_with_rod(this, position, orientation, res, err)
       import dp
       import rod_interaction
       class(rod_interaction), intent(inout) :: this
       real(dp), intent(in) :: position(3), orientation(3)
       real(dp), intent(out) :: res
       integer, intent(out) :: err
     end subroutine interact_with_rod
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

  elemental subroutine rod_energy(this, ia, res, err)
    class(rod), intent(in) :: this
    class(rod_interaction), intent(inout) :: ia
    real(REAL64), intent(out) :: res
    integer, intent(out) :: err
    call ia%value(this%position, this%orientation, res, err)
  end subroutine

end module m_rod



module m_particle_interaction_lists
  use m_particle, only: particle_interaction
  type particle_interaction_ptr
     class(particle_interaction), pointer :: ptr
  end type particle_interaction_ptr
  
#define TYPEPARAM type(particle_interaction_ptr)
#include "list-inc-def.inc"
contains
#include "list-inc-proc.inc"
#undef TYPEPARAM
end module m_particle_interaction_lists



module m_particlegroup
  use nrtype
  use m_particle
  use json_module
  use m_particle_interaction_lists, only: &
       particle_ia_list => list, &
       particle_interaction_ptr
  implicit none

  type particlegroup
     character(kind=CK, len=:), allocatable :: name
     real(dp), allocatable :: positions(:, :)
     type(particle_ia_list) :: particle_interactions
   contains
     procedure :: add_interaction => add_gen_ia
     !   final :: particle_finalize
     procedure :: to_json => particlegroup_to_json
     procedure :: type_str => particlegroup_type_str
     !procedure(reduce), deferred :: reduce_single
     procedure :: reduce
  end type particlegroup
  
  interface particlegroup
     module procedure empty_group, particlegroup_from_json
  end interface particlegroup

  !interface
  !   subroutine reduce(this, potential, res, err)
  !     import particlegroup, particle_interaction, REAL64
  !     class(particlegroup), intent(in) :: this
  !     class(particle_interaction), intent(in) :: potential
  !     real(REAL64), intent(out) :: res
  !     integer, intent(out) :: err
  !   end subroutine reduce
  !end interface
  
contains 
  
  function empty_group(name)
    type(particlegroup) :: empty_group
    character(len=:), allocatable, intent(in) :: name
    empty_group%name = name(1:80)
  end function empty_group

  function particlegroup_from_json(json_val) result(pg)
    type(particlegroup) :: pg
    type(json_value), pointer, intent(in) :: json_val
    type(json_value), pointer :: child => null()
    logical :: found
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
    character(kind=CK, len=:), allocatable :: str
    call json_create_object(json_val, '')
    call json_add(json_val, 'name', this%name)
    call this%type_str(str)
    call json_add(json_val, 'type', str)
    call add_3d_points_json(json_val, 'positions', this%positions)
  end subroutine particlegroup_to_json

  subroutine particlegroup_type_str(this, str)
    class(particlegroup), intent(in) :: this
    character(kind=CK, len=:), allocatable :: str
    str = 'particlegroup'
  end subroutine particlegroup_type_str
  
  subroutine add_3d_points_json(json_val, name, points)
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=*), intent(in) :: name
    real(REAL64), intent(in) :: points(:, :)
    type(json_value), pointer :: temp_array => null()
    integer :: i
    call json_create_array(temp_array, name)
    do i = 1, size(points(1, :)) 
       call json_add(temp_array, '', points(:, i))
    end do
    call json_add(json_val, temp_array)
  end subroutine add_3d_points_json

  subroutine add_gen_ia(this, ia)
    class(particlegroup), intent(inout) :: this
    class(particle_interaction), intent(in), target :: ia
    call this%particle_interactions%add(particle_interaction_ptr(ia))
  end subroutine add_gen_ia
  
  subroutine energy(this, aparticle, res, err)
    class(particlegroup), intent(inout) :: this
    class(particle), intent(in) :: aparticle
    real(dp), intent(out) :: res
    integer, intent(out) :: err
    type(particle_interaction_ptr), pointer :: ia
    real(dp) :: cur
    res = 0.
    call this%particle_interactions%iter_restart()
    do while(.true.)
       ia => this%particle_interactions%iter_next()
       if (.not. associated(ia)) exit
       call aparticle%energy(ia%ptr, cur, err)
       if (err /= 0) exit
       res = res + cur
    end do
  end subroutine energy
  
  subroutine reduce(this, potential, res, err)
    class(particlegroup), intent(in) :: this
    class(particle_interaction), intent(in) :: potential
    real(REAL64), intent(out) :: res
    logical, intent(out) :: err
  end subroutine reduce

  subroutine particle_finalize(this)
    type(particlegroup), intent(inout) :: this
    call this%particle_interactions%finalize()
  end subroutine particle_finalize
  
end module m_particlegroup


!!example of parametric list usage
!#define STRING32 32
module m_rod_interaction_lists
  use m_rod
  type rod_interaction_ptr
     class(rod_interaction), pointer :: ptr
  end type rod_interaction_ptr
  
#define TYPEPARAM type(rod_interaction_ptr)
#include "list-inc-def.inc"
contains
#include "list-inc-proc.inc"
#undef TYPEPARAM
end module


module m_rodgroup
  
  use m_rod
  use m_particlegroup
  use m_rod_interaction_lists, only: rod_ia_list => list
  implicit none

  !type, extends(particlegroup) :: rodgroup
  type rodgroup
     type(rod_ia_list) :: interactions
     real(dp), allocatable :: orientations(:, :)
   contains
!     procedure :: add_interaction => add_rod_ia
     procedure :: to_json => rodgroup_to_json
     procedure :: type_str => rodgroup_type_str
     procedure :: reduce => reduce_rods
  end type rodgroup

  interface rodgroup
     module procedure rodgroup_from_json
  end interface rodgroup

contains

  function rodgroup_from_json(json_val) result(rg)
    type(rodgroup) :: rg
    type(json_value), pointer, intent(in) :: json_val
    type(json_value), pointer :: child
    logical :: found
    !rg%particlegroup = particlegroup(json_val)
    call json_get(json_val, 'orientations', child, found)
    call get_3d_points_json(child, rg%orientations)
  end function rodgroup_from_json

  subroutine rodgroup_to_json(this, json_val)
    class(rodgroup), intent(in) :: this
    type(json_value), pointer :: json_val
    character(kind=CK, len=:), allocatable :: str
    logical :: found
    !call this%particlegroup%to_json(json_val)
    call add_3d_points_json(json_val, 'orientations', this%orientations)
    call this%type_str(str)
    call json_update(json_val, 'type', str, found)
  end subroutine rodgroup_to_json

  subroutine rodgroup_type_str(this, str)
    class(rodgroup), intent(in) :: this
    character(kind=CK, len=:), allocatable :: str
    str = 'rodgroup'
  end subroutine rodgroup_type_str
  
  subroutine add_rod_ia(this, ia)
    class(rodgroup), intent(inout) :: this
    class(rod_interaction), intent(in) :: ia
  end subroutine add_rod_ia

  pure subroutine reduce_rods(this, potential, res, err)
    class(rodgroup), intent(in) :: this
    class(rod_interaction), intent(in) :: potential
    real(REAL64), intent(out) :: res
    logical, intent(out) :: err
  end subroutine reduce_rods
  
end module m_rodgroup



!! define a particlegroup list
module m_particlegroup_list
  use m_particlegroup
#define TYPEPARAM type(particlegroup)
#include "list-inc-def.inc"
contains
#include "list-inc-proc.inc"
#undef TYPEPARAM
end module m_particlegroup_list


!! define a rodgroup list
module m_rodgroup_list
  use m_rodgroup
#define TYPEPARAM type(rodgroup)
#include "list-inc-def.inc"
contains
#include "list-inc-proc.inc"
#undef TYPEPARAM
end module m_rodgroup_list


module m_particleserver

  use m_particlegroup
  use m_rodgroup
  use m_rodgroup_list, only: &
       rodgroup_list => list, &
       add_rodgroup => list_add, &
       rodgroup_list_restart => list__iter_restart, &
       next_rodgroup => list_iter_next
  use m_particlegroup_list, only: &
       particlegroup_list => list, &
       add_particlegroup => list_add
  implicit none
  
  type particleserver_json
     type(rodgroup_list) :: rodgroups
     type(particlegroup_list) :: particlegroups
   contains
     !  procedure :: get_rodgroup
     procedure :: build
     procedure :: serialize
  end type particleserver_json

  interface particleserver_json
     module procedure particleserver_from_json
  end interface particleserver_json
  
contains

  function particleserver_from_json(filename) result(ps)
    character(kind=CK, len=*), intent(in) :: filename
    type(particleserver_json) :: ps
    call ps%build(filename)
  end function particleserver_from_json
  
  subroutine build(this, filename)
    class(particleserver_json), intent(inout) :: this
    character(len=*), intent(in) :: filename
    type(json_file) :: json
    logical :: found
    type(json_value), pointer :: p, child
    integer(INT32) :: i
    call json_initialize()

    ! read particlegroups from file        
    call json%load_file(filename = filename)!'groups.json')
    call json%get('groups', p, found)
    !allocate(groups(json_count(p)))
    do i = 1, json_count(p)
       !write(output_unit, *) 'Get particle from json:'
       call json_get_child(p, i, child)
       !write(output_unit, *) 'Read particle from json:'
       ! build particlegrops
       call particle_factory(this, child)
    end do
    
    ! clean up
    call json%destroy()
    if (json_failed()) stop 1
    
  end subroutine build

  subroutine particle_factory(this, json_val)
    class(particleserver_json), intent(inout) :: this
    type(json_value), pointer, intent(in) :: json_val
    character(kind=CK, len=:), allocatable :: type_str
    logical :: found
    call json_get(json_val, 'type', type_str, found)
    if (.not. found) then
       write(error_unit, *) 'Error: Could not read particle type!'
       stop 1
    end if
    type_str = adjustl(type_str)
    if (trim(type_str) == "particlegroup") then
       call this%particlegroups%add(particlegroup(json_val))
    else if (trim(type_str) == "rodgroup") then
       call this%rodgroups%add(rodgroup(json_val))
    end if
  end subroutine particle_factory

  
  subroutine serialize(this, filename)
    class(particleserver_json), intent(inout) :: this
    character(len=*), intent(in) :: filename
    class(particlegroup), pointer :: res
    integer :: n = 0
    type(json_value), pointer :: json_val, group_json
    character(kind=CK, len=:), allocatable :: str
    nullify(res, json_val, group_json)
    call json_create_array(json_val, 'groups')

    call this%rodgroups%iter_restart()
    do while(.true.)
       res => this%rodgroups%iter_next()
       if (.not. associated(res)) then
          exit
       end if
       n = n + 1
       !select type (res)
       !class is (particlegroup)
       call res%to_json(group_json)
       call json_add(json_val, group_json)
       !end select
    end do

    call this%particlegroups%iter_restart()
    do while(.true.)
       res => this%particlegroups%iter_next()
       if (.not. associated(res)) then
          exit
       end if
       n = n + 1
       !select type (res)
       !class is (particlegroup)
       call res%to_json(group_json)
       call json_add(json_val, group_json)
       !end select
    end do
    call json_print_to_string(json_val, str)
    write(*, *) str
  end subroutine serialize
  
  function get_rodgroup(this, name, found) result(rg)
    class(particleserver_json), intent(inout) :: this
    character(kind=CK, len=*), intent(in) :: name
    logical, intent(out) :: found
    type(rodgroup), pointer :: rg
    nullify(rg)
    call this%rodgroups%iter_restart()
    do while(.true.)
       rg => this%rodgroups%iter_next()
       if (.not. associated(rg)) then
          found = .false.
          exit
       end if
       if (rg%name == 'name') return
    end do
  end function get_rodgroup
  
end module m_particleserver




program test_composition
  use json_module
  use m_particleserver
  implicit none
  !class(particle), allocatable :: groups(:)
  !type(list) :: groups
  type(particleserver_json) :: particleserver
  ! Read particlegroups from file
  !call build('groups.json', groups)
  particleserver = particleserver_json('groups.json')
  
  ! read interactions from file
  !interactionserver = interactionserver_json('interactions.json')

  ! attach interactions to particlegroups
  !call interactionserver%connect_groups(particleserver)  
  ! simulate

  ! write results to file
  call particleserver%serialize('output.json')
  !call interactionserver%serialize('output.json')
end program test_composition

