!> Contains the point type and related procedures.
module m_point
  use m_particle, only: point
  use num_kind, only: dp
  use particle_mover, only: transmove, rotate
  use json_module, only: json_value, json_add, json_create_array, CK
  use mt_stream, only: rng=>genrand_double1_s, rngstate=>mt_state
  use m_json_wrapper, only: get_parameter
  use utils, only: fmt_char_dp
  implicit none

  !> Point-like particle.
  type, extends(point) :: point2
   contains
     !> Deserializes the coordinates from a string.
     procedure :: from_str => point_from_str
     !> Downcasts a particle to a point and assigns
     procedure :: downcast_assign => point_downcast_assign
     !> Assignment operator implementation
     procedure :: point_assign
     generic :: assignment(=) => point_assign
     !> Equality operator implementation.
     procedure :: point_equals
     generic :: operator(==) => point_equals
     !> Returns the type of this particle as a string.
     procedure, nopass :: typestr => point_typestr
     !> Randomly translates the point to a new position.
     procedure :: move => point_move
     !> Writes point coordinates to a JSON value.
     procedure :: coordinates_to_json => point_coordinates_to_json
     !> Deserializes point coordinates from json.
     procedure :: from_json => point_from_json
     !> Returns the names of the point coordinates.
     procedure, nopass :: description => point_description
     !> Serializes the point to the given Fortran output unit.
     procedure :: to_unit => point_to_unit
  end type point2


contains

  
  !> Deserializes @p this point from @p str. ios /= 0 if an error
  !! occurs.
  subroutine point_from_str(this, str, ios)
    class(point), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z
  end subroutine point_from_str


  !> Returns the type of the point in @p str.
  subroutine point_typestr(str)
    character(len=:), allocatable, intent(out) :: str
    str = "point"
  end subroutine point_typestr


  !> Returns the names of the point coordinates in @p descr.
  subroutine point_description(descr)
    character(kind=CK, len=3), allocatable, intent(inout) :: descr(:)
    descr = ["x  ", "y  ", "z  "]
  end subroutine point_description


  !> Downcasts @p a_particle to a point and assigns it to @p this.
  !! err = 3 if an error occurs.
  pure subroutine point_downcast_assign(this, a_particle, err)
    class(point), intent(inout) :: this
    class(point), intent(in) :: a_particle
    integer, intent(out), optional :: err
    select type (a_particle)
    type is (point)
       this = a_particle
       if (present(err)) err = 0
    class default
       if (present(err)) err = 3
    end select
  end subroutine point_downcast_assign


  !> Assigns @p another to @p this.
  pure subroutine point_assign(this, another)
    class(point), intent(inout) :: this
    type(point), intent(in) :: another
    this%x = another%x
    this%y = another%y
    this%z = another%z
  end subroutine point_assign


  !> Returns true if coordinates of @p another are equal to @p this.
  elemental function point_equals(this, another) result(res)
    class(point), intent(in) :: this
    type(point), intent(in) :: another
    logical :: res
    res = (this%x == another%x) .and. (this%y == another%y) .and. &
         (this%z == another%z)
  end function point_equals


  !> Translates @p this randomly to a new position. @p genstate is the
  !! random number generator used.
  pure subroutine point_move(this, genstate)    
    class(point), intent(inout) :: this
    type(rngstate), intent(inout) :: genstate
    real(dp) :: xn, yn, zn
    call transmove(this%x, this%y, this%z, xn, yn, zn, genstate)
    this%x = xn
    this%y = yn
    this%z = zn
  end subroutine point_move


  !> Serializes the coordinates of @p this point to JSON value
  !! @p json_val.
  subroutine point_coordinates_to_json(this, json_val)
    class(point), intent(in) :: this
    type(json_value), pointer :: json_val
    call json_create_array(json_val, '')
    call json_add(json_val, '', this%x)
    call json_add(json_val, '', this%y)
    call json_add(json_val, '', this%z)
  end subroutine point_coordinates_to_json


  !> Deserializes @p this point from JSON value @p json_val.
  subroutine point_from_json(this, json_val)
    class(point), intent(inout) :: this
    type(json_value), pointer, intent(in) :: json_val
    call get_parameter(json_val, '[1]', this%x)
    call get_parameter(json_val, '[2]', this%y)
    call get_parameter(json_val, '[3]', this%z)
  end subroutine point_from_json


  !> Writes @p this point to @p unit.
  subroutine point_to_unit(this, unit)
    class(point), intent(in) :: this
    integer, intent(in) :: unit
    write(unit, '(A,3(' // fmt_char_dp() // ',1X))', advance='no') &
         'point', this%x, this%y, this%z
  end subroutine point_to_unit

end module m_point
