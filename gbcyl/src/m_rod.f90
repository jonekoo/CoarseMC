module m_rod
  use iso_fortran_env, only: output_unit
  use particle, only: particledat
  use num_kind, only: dp
  use utils, only: fmt_char_dp
  use particle_mover, only: transmove, rotate
  include 'rng.inc'
  implicit none
  
  type, extends(particledat) :: rod
   contains
     procedure :: rod_from_str
     procedure :: downcast_assign => rod_downcast_assign
     procedure :: rod_assign
     generic :: assignment(=) => rod_assign
     procedure :: rod_equals
     generic :: operator(==) => rod_equals
     procedure, nopass :: typestr => rod_typestr
     procedure :: to_stdout => rod_to_stdout
     procedure :: move => rod_move
  end type rod

contains

  subroutine rod_from_str(this, str, ios)
    class(rod), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z, this%ux, &
         this%uy, this%uz
  end subroutine rod_from_str

  subroutine rod_typestr(str)
    character(len=:), allocatable, intent(out) :: str
    str = "rod"
  end subroutine rod_typestr

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

  elemental function rod_equals(this, another) result(res)
    class(rod), intent(in) :: this
    type(rod), intent(in) :: another
    logical :: res
    res = (this%x == another%x) .and. (this%y == another%y) .and. &
         (this%z == another%z) .and. (this%ux == another%ux) .and. &
         (this%uy == another%uy) .and. (this%uz == another%uz)
  end function rod_equals
  
  subroutine rod_to_stdout(this)
    class(rod), intent(in) :: this
    write(output_unit, '(6(' // fmt_char_dp() // ',1X))', advance='no') &
         this%x, this%y, this%z, this%ux, this%uy, this%uz 
  end subroutine rod_to_stdout

  pure subroutine rod_move(this, genstate)    
    class(rod), intent(inout) :: this
    type(rngstate), intent(inout) :: genstate
    real(dp) :: xn, yn, zn, uxn, uyn, uzn
    call transmove(this%x, this%y, this%z, xn, yn, zn, genstate)
    this%x = xn
    this%y = yn
    this%z = zn
    call rotate(this%ux, this%uy, this%uz, uxn, uyn, uzn, genstate)
    this%ux = uxn
    this%uy = uyn
    this%uz = uzn
  end subroutine rod_move

end module m_rod
