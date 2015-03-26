module m_base
type base
contains
procedure :: print_type
end type

contains

subroutine print_type(this)
  class(base), intent(in) :: this
  write(*, *) 'with base'
end subroutine

end module


module m_child
use m_base
type, extends(base) :: child
contains
procedure :: print_type => print_child
end type

contains

subroutine print_child(this)
  class(child), intent(in) :: this
  write(*, *) 'with child'
end subroutine

end module


module m_operation
use m_base
use m_child

type operation
  contains
    procedure :: for_base
    generic :: for_who => for_base
end type

contains

subroutine for_base(this, abase)
  class(operation), intent(in) :: this
  class(base), intent(in) :: abase
  write(*, *) 'for_base called'
  call abase%print_type()
end subroutine

end module






program test
use m_base
use m_child
use m_operation
class(base), allocatable :: abase, achild
class(operation), allocatable :: op
allocate(base::abase)
allocate(child::achild)
allocate(operation::op)
call op%for_who(abase)
call op%for_who(achild)
end program
