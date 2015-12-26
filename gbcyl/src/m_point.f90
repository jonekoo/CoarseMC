module m_point
  use particle, only: particledat
  use num_kind, only: dp
  use particle_mover, only: transmove, rotate
  include 'rng.inc'
  implicit none

  type, extends(particledat) :: point
   contains
     procedure :: from_str => point_from_str
     procedure :: downcast_assign => point_downcast_assign
     procedure :: point_assign
     generic :: assignment(=) => point_assign
     procedure :: point_equals
     generic :: operator(==) => point_equals
     procedure, nopass :: typestr => point_typestr
     procedure :: move => point_move
  end type point

contains
  
  subroutine point_from_str(this, str, ios)
    class(point), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z
  end subroutine point_from_str
  
  subroutine point_typestr(str)
    character(len=:), allocatable, intent(out) :: str
    str = "point"
  end subroutine point_typestr
  
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

  elemental function point_equals(this, another) result(res)
    class(point), intent(in) :: this
    type(point), intent(in) :: another
    logical :: res
    res = (this%x == another%x) .and. (this%y == another%y) .and. &
         (this%z == another%z)
  end function point_equals
  
  pure subroutine point_move(this, genstate)    
    class(point), intent(inout) :: this
    type(rngstate), intent(inout) :: genstate
    real(dp) :: xn, yn, zn, uxn, uyn, uzn
    call transmove(this%x, this%y, this%z, xn, yn, zn, genstate)
    this%x = xn
    this%y = yn
    this%z = zn
  end subroutine point_move

end module m_point
