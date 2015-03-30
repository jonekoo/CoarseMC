module m_particle
  use nrtype
  use gen_lists

type particle
   type(list) :: interactions
   real(dp) :: position(3) = [0., 0., 0.]
   contains
   procedure :: add_interaction => add_gen_ia
!   final :: particle_finalize
end type particle

contains

subroutine add_gen_ia(this, ia)
  class(particle), intent(inout) :: this
  class(*), intent(in) :: ia
  call this%interactions%add(ia)
end subroutine

subroutine particle_finalize(this)
  type(particle), intent(inout) :: this
  call this%interactions%finalize()
end subroutine

end module



module m_rod
  use m_particle
  use iso_fortran_env

type, extends(particle) :: rod
  real(dp) :: orientation(3) = [0., 0., 1.]
  contains
  procedure :: add_interaction => add_rod_ia
end type

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
end subroutine     

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
end subroutine

end module



module m_particlegroup
  use nrtype
  use gen_lists
  use m_particle

type, extends(particle) :: particlegroup
  character(len=:), allocatable :: name
  real(dp), allocatable :: positions(:, :)
end type

interface particlegroup
   module procedure empty_group
end interface

contains 

function empty_group(name)
  type(particlegroup) :: empty_group
  character(len=:), allocatable, intent(in) :: name
  empty_group%name = name(1:80)
end function

end module

module m_rodgroup
  use m_rod

type, extends(rod) :: rodgroup
  real(dp), allocatable :: orientations(:, :)
  real(dp), allocatable :: positions(:, :)
end type

end module


module m_builder

  use m_particlegroup
  use iso_fortran_env
  use json_module


contains 


subroutine build(filename, groups)
  character(len=*), intent(in) :: filename
  class(particle), allocatable, intent(out) :: groups(:)
  type(json_file) :: json
  logical :: found
  !integer :: i,j,k
  type(json_value), pointer :: p, child, pos
  integer(INT32) :: i
  character(kind=CK, len=1) :: ic
  character(kind=CK, len=:), allocatable :: name
  character(kind=CK, len=:), allocatable :: type
  integer(INT32) :: group_size
  real(REAL64), allocatable :: positions(:, :)
  real(REAL64), allocatable :: temp(:)
  call json_initialize()

  ! read particlegroups from file        
  call json%load_file(filename = filename)!'groups.json')
  call json%get('groups', p, found)
  allocate(groups(json_count(p)))
  do i = 1, json_count(p) 
     call json_get_child(p, i, child)

     call json_get(child, 'name', name, found)
     if (.not. found) then
        stop 1
     end if


     call json_get(child, 'size', group_size, found)
     if (.not. found) then
        stop 1
     end if

     call json_get(child, 'positions', pos, found)
     if (.not. found) then
        stop 1
     end if
     allocate(positions(3, json_count(pos)))
     do j = 1, json_count(pos)
        call json_get(pos, j, temp)
     end do

     call json_get(child, 'type', type, found)
     if (.not. found) then
        stop 1
     end if
     write(output_unit, *) name, group_size, type
     !if (type == 'rodgroup') then
   
        !groups(1) = rodgroup(

  end do
  
  ! clean up
  call json%destroy()
  if (json_failed()) stop 1



  ! build particlegrops
end subroutine build

end module

program test_composition
  use m_builder
  ! read interactions from file
  ! attach interactions to particlegroups
  class(particle), allocatable :: groups(:)
  call build('groups.json', groups)
  ! simulate
end program

