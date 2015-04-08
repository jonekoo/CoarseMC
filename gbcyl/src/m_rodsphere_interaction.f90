module m_rodsphere_interaction
use m_particlegroup
use m_rodgroup

type, abstract :: rodsphere_potential
 contains
   procedure(rodsphere_value), deferred :: value
end type rodsphere_potential

interface
   subroutine rodsphere_value(this, u, r, res, err)
     import rodsphere_potential, REAL64
     class(rodsphere_potential), intent(in) :: this
     real(REAL64), intent(in) :: u(3)
     real(REAL64), intent(in) :: r(3)
     real(REAL64), intent(out) :: res
     integer, intent(out) :: err
   end subroutine rodsphere_value
end interface


!! This could actually be a general particle-rod interaction with a
!! variable component to compute the potential.
type rodsphere_interaction
   character(kind=CK, len=:), allocatable :: rodgroup_name
   character(kind=CK, len=:), allocatable :: particlegroup_name
   class(rodsphere_potential), pointer :: potential => null()
end type rodsphere_interaction

!! These could be implemented as wrappers around the general particle-rod
!! interaction.
type, extends(particle_interaction) ::  rodsphere_bind_particle
  type(particle) :: p
  class(rodgroup), pointer :: rg => null()
  class(rodsphere_potential), pointer :: potential
  contains 
    procedure :: value => value_particle
end type

type, extends(rod_interaction) :: rodsphere_bind_rod
  type(rod) :: r
  class(particlegroup), pointer :: pg => null()
  class(rodsphere_potential), pointer :: potential
  contains
    procedure :: value => value_rod
end type

contains

  subroutine connect(this, groups)
    class(rodsphere_interaction), intent(inout) :: this
    class(particlegroup), intent(in), target :: groups(:)
    class(rodgroup), pointer :: rg
    class(particlegroup), pointer :: pg
    nullify(rg, pg)
    do i = 1, size(groups)
       if (groups(i)%name == this%rodgroup_name) then
          select type (groups(i))
          class is (rodgroup)
             rg => groups(i)
          end select
       end if
       if (.not. associated(rg)) then
          write(error_unit, *) "Error: rodpshere_interaction: rodgroup with " &
               // "name " // trim(this%rodgroup_name) // " not found. Stopping."
          stop
       end if
       if (groups(i)%name == this%particlegroup_name) then
          select type (groups(i))
          class is (particlegroup)
             pg => groups(i)
          end select
       end if
       if (.not. associated(pg)) then
          write(error_unit, *) "Error: rodpshere_interaction: particlegroup " &
               // "with name " // trim(this%particlegroup_name) // &
               " not found. Stopping."
          stop
       end if
    end do
    call rg%add_interaction(rodsphere_bind_rod(pg, this%potential))
    call pg%add_interaction(rodsphere_bind_particle(rg, this%potential))
end subroutine connect


function rodsphere_bind_rod_pg(pg, potential) result(this)
  class(particlegroup), pointer, intent(in) :: pg
  class(rodsphere_potential), pointer, intent(in) :: potential
  type(rodsphere_bind_rod) :: this
  nullify(this%pg)
  if (associated(pg)) then
     this%pg => pg
  else
     write(error_unit, *) 'Error: rodsphere_bind_rod_pg: pointer pg is null.'
  end if
  this%potential => potential
end function rodsphere_bind_rod_pg


function rodsphere_bind_particle_rg(rg, potential) result(this)
  class(rodgroup), pointer, intent(in) :: rg
  class(rodsphere_potential), pointer, intent(in) :: potential
  type(rodsphere_bind_particle) :: this
  nullify(this%rg)
  if (associated(rg)) then
     this%rg => rg
  else
     write(error_unit, *) 'Error: rodsphere_bind_particle_rg(rg, potential)'
  end if
  this%potential => potential
end function rodsphere_bind_particle_rg

pure subroutine value_particle(this, position, res, err)
  class(rodsphere_bind_particle), intent(inout) :: this
  real(REAL64), intent(in) :: position(3)
  real(REAL64), intent(out) :: res
  integer, intent(out) :: err
  this%p = p
  call this%rg%reduce(self, res, err)
end subroutine value_particle

pure subroutine value_rod(this, position, orientation, res, err)
  class(rodsphere_bind_rod), intent(inout) :: this
  real(REAL64), intent(in) :: position(3), orientation(3)
  real(REAL64), intent(out) :: res
  integer, intent(out) :: err
  this%r = rod(position, orientation)
  call this%pg%reduce(this, res, err)
end subroutine value_rod

end module m_rodsphere_interaction
