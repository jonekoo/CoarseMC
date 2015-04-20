module m_types
  implicit none
  
  type foo
     integer :: i = 1
   contains
     procedure :: modify => modify_i
     procedure :: str => foo_str
  end type foo

  type, extends(foo) :: bar
     integer :: j = 2
   contains
     procedure :: modify => modify_j
     procedure :: str => bar_str
     procedure :: bar_assign_foo
     generic :: assignment(=) => bar_assign_foo
  end type bar


contains

  subroutine modify_i(this)
    class(foo), intent(inout) :: this
    this%i = this%i + 1
  end subroutine modify_i

  function foo_str(this)
    class(foo), intent(in) :: this
    character(len=:), allocatable :: foo_str
    allocate(character(5) :: foo_str)
    write(foo_str, fmt='(A, 1X, I1)') 'foo', this%i
  end function foo_str
  
  subroutine modify_j(this)
    class(bar), intent(inout) :: this
    this%j = this%j + 1
  end subroutine modify_j

  function bar_str(this)
    class(bar), intent(in) :: this
    character(len=:), allocatable :: bar_str
    allocate(character(7) :: bar_str)
    write(bar_str, fmt='(A, 1X, I1, 1X, I1)') 'bar', this%i, this%j
  end function bar_str

  subroutine bar_assign_foo(this, src)
    class(foo), intent(in) :: src
    class(bar), intent(out) :: this
    this%i = src%i
    select type (src)
    class is (bar)
       this%j = src%j
    end select
  end subroutine bar_assign_foo
  
end module m_types



program polymorphic_copy
  use iso_fortran_env
  use m_types
  type(bar) :: modifiable

  call modify()
  
contains
  
  subroutine modify()
    class(foo), allocatable :: temp_copy
    allocate(temp_copy, source=modifiable)
    call temp_copy%modify()
    write(output_unit, *) 'After modify:'
    write(output_unit, *) 'modifiable = ', modifiable%str()
    write(output_unit, *) 'temp_copy = ', temp_copy%str()
    modifiable = temp_copy
    write(output_unit, *) 'After copying back:'
    write(output_unit, *) 'modifiable = ', modifiable%str()
    write(output_unit, *) 'temp_copy = ', temp_copy%str()
  end subroutine
  
end program polymorphic_copy
