module m_point
  use particle, only: particledat
  use num_kind, only: dp
  implicit none

  type, extends(particledat) :: point
   contains
     procedure :: from_str => point_from_str
  end type point

contains
  
  subroutine point_from_str(this, str, ios)
    class(point), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z
  end subroutine point_from_str
  
end module m_point
