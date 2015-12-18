module m_rod
  use particle, only: particledat
  use num_kind, only: dp
  implicit none
  
  type, extends(particledat) :: rod
   contains
     procedure :: rod_from_str
     procedure :: downcast_assign => rod_downcast_assign
     procedure :: rod_assign
     generic :: assignment(=) => rod_assign
  end type rod

contains

  subroutine rod_from_str(this, str, ios)
    class(rod), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z, this%ux, &
         this%uy, this%uz
  end subroutine rod_from_str
  
  pure subroutine rod_downcast_assign(this, a_particle, err)
    class(rod), intent(inout) :: this
    class(particledat), intent(in) :: a_particle
    integer, intent(out), optional :: err
    select type (a_particle)
    type is (rod)
       this = a_particle
    class default
       if (present(err)) err = 3
    end select
  end subroutine rod_downcast_assign

  pure subroutine rod_assign(this, another)
    class(rod), intent(inout) :: this
    type(rod), intent(in) :: another
    this%x = another%x
    this%y = another%y
    this%z = another%z
    this%ux = another%ux
    this%uy = another%uy
    this%uz = another%uz
  end subroutine rod_assign

end module m_rod
