module m_rod
  use particle, only: particledat
  use num_kind, only: dp
  implicit none
  
  type, extends(particledat) :: rod
   contains
     procedure :: rod_from_str
  end type rod

contains

  subroutine rod_from_str(this, str, ios)
    class(rod), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z, this%ux, &
         this%uy, this%uz
  end subroutine rod_from_str
  
end module m_rod
