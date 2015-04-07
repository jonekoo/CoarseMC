module m_gblj_interaction
use m_particlegroup
use m_rodgroup

!! This could actually be a general particle-rod interaction with a
!! variable component to compute the potential.
type gblj_interaction
   character(kind=CK, len=:), allocatable :: rodgroup_name = ""
   character(kind=CK, len=:), allocatable :: particlegroup_name = ""
   class(rodgroup), pointer :: rg => null()
   class(particlegroup), pointer :: pg => null()
end type gblj_interaction

!! These could be implemented as wrappers around the general particle-rod
!! interaction.
type, extends(particle_interaction) ::  gblj_bind_particle
  type(particle) :: p
  class(rodgroup), pointer :: rg => null()
  contains 
    procedure :: value => value_particle
end type

type, extends(rod_interaction) :: gblj_bind_rod
  type(rod) :: r
  class(particlegroup), pointer :: pg => null()
  contains
    procedure :: value => value_rod
end type

contains

subroutine connect(this, groups)
  class(*), intent(in), target :: groups(:)
  do i = 1, size(groups)
    if (groups(i)%name == this%rodgroup_name) then
        select type (groups(i))
        class is (rodgroup)
          this%rg => groups(i)
       end select
    end if
    if (.not. associated(this%rg)) then
       write(error_unit, *) "Error: gblj_interaction: rodgroup with name " //&
            trim(this%rodgroup_name) // " not found. Stopping."
       stop
    end if
    if (groups(i)%name == this%particlegroup_name) then
       select type (groups(i))
       class is (particlegroup)
         this%pg => groups(i)
      end select
    end if
    if (.not. associated(this%pg)) then
      write(error_unit, *) "Error: gblj_interaction: particlegroup with name "&
           // trim(this%particlegroup_name) // " not found. Stopping."
      stop
    end if
  end do
  this%rg%add_interaction(gblj_bind_rod(this%pg))
  this%pg%add_interaction(gblj_bind_particle(this%rg))
end subroutine

end module
