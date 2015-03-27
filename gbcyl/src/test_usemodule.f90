module points2d

type point2d
   real :: x, y
contains
  procedure :: print_name
end type

contains

subroutine print_name(this)
  class(point2d), intent(in) :: this
  write(*, *) 'original point!'
end subroutine print_name

end module points2d




module new_points2d
use points2d, point2d_original => point2d, print_name_original => print_name

type, extends(point2d_original) :: point2d
contains
  procedure :: print_name 
end type

contains 

subroutine print_name(this)
  class(point2d), intent(in) :: this
  call this%point2d_original%print_name()
  write(*, *) 'new point!'
end subroutine print_name

end module new_points2d

module point2d_functionality
use new_points2d
!use points2d
end module


program test_usemodule
use point2d_functionality
class(point2d), allocatable :: p
allocate(point2d::p)
call p%print_name()
end program
