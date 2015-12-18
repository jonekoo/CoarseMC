module m_point
  use particle, only: particledat
  use num_kind, only: dp
  implicit none

  type, extends(particledat) :: point
   contains
     procedure :: from_str => point_from_str
     procedure :: downcast_assign => point_downcast_assign
     procedure :: point_assign
     generic :: assignment(=) => point_assign
  end type point

contains
  
  subroutine point_from_str(this, str, ios)
    class(point), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z
  end subroutine point_from_str
  
  pure subroutine point_downcast_assign(this, a_particle, err)
    class(point), intent(inout) :: this
    class(particledat), intent(in) :: a_particle
    integer, intent(out), optional :: err
    select type (a_particle)
    type is (point)
       this = a_particle
    class default
       if (present(err)) err = 3
    end select
  end subroutine point_downcast_assign

  pure subroutine point_assign(this, another)
    class(point), intent(inout) :: this
    type(point), intent(in) :: another
    this%x = another%x
    this%y = another%y
    this%z = another%z
  end subroutine point_assign

end module m_point
