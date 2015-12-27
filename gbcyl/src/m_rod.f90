module m_rod
  use iso_fortran_env, only: output_unit
  use particle, only: particledat
  use num_kind, only: dp
  use utils, only: fmt_char_dp
  use particle_mover, only: transmove, rotate
  include 'rng.inc'
  use json_module, only: json_create_array, CK, json_value, json_add
  use m_json_wrapper, only: get_parameter
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
     procedure :: coordinates_to_json => rod_coordinates_to_json
     procedure :: from_json => rod_from_json
     procedure, nopass :: description => rod_description
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

  subroutine rod_description(descr)
    character(kind=CK, len=3), allocatable, intent(inout) :: descr(:)
    descr = ["x  ", "y  ", "z  ", "ux ", "uy ", "uz "]
  end subroutine rod_description

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

  subroutine rod_coordinates_to_json(this, json_val)
    class(rod), intent(in) :: this
    type(json_value), pointer :: json_val
    call json_create_array(json_val, '')
    call json_add(json_val, '', this%x)
    call json_add(json_val, '', this%y)
    call json_add(json_val, '', this%z)
    call json_add(json_val, '', this%ux)
    call json_add(json_val, '', this%uy)
    call json_add(json_val, '', this%uz)
  end subroutine rod_coordinates_to_json
  
  subroutine rod_from_json(this, json_val)
    class(rod), intent(inout) :: this
    type(json_value), pointer, intent(in) :: json_val
    call get_parameter(json_val, '[1]', this%x)
    call get_parameter(json_val, '[2]', this%y)
    call get_parameter(json_val, '[3]', this%z)
    call get_parameter(json_val, '[4]', this%ux, error_lb=-1._dp, &
         error_ub=1._dp)
    call get_parameter(json_val, '[5]', this%uy, error_lb=-1._dp, &
         error_ub=1._dp)
    call get_parameter(json_val, '[6]', this%uz, error_lb=-1._dp, &
         error_ub=1._dp)
  end subroutine rod_from_json

end module m_rod
