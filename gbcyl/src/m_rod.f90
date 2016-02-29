!> Implements a rod-like particle.
module m_rod
  use iso_fortran_env, only: output_unit
  use m_particle, only: particle
  use m_point, only: point
  use num_kind, only: dp
  use utils, only: fmt_char_dp
  use particle_mover, only: transmove, rotate
  use json_module, only: json_create_array, CK, json_value, json_add
  use m_json_wrapper, only: get_parameter
  use mt_stream, only: rng=>genrand_double1_s, rngstate=>mt_state
  implicit none

  !> A rod-like particle.
  type, extends(point) :: rod
     !> The components of the unit vector of orientation.
     real(dp) :: ux = 0._dp
     real(dp) :: uy = 0._dp
     real(dp) :: uz = 1._dp
   contains
     !> Returns the type of the rod as a string.
     procedure, nopass :: typestr => rod_typestr
     !> Returns the description of the coordinates for this rod as an
     !! array of strings.
     procedure, nopass :: description => rod_description

     !> Assignment operator implementation.
     procedure :: rod_assign
     generic :: assignment(=) => rod_assign
     !> Downcasts and assigns a particle to a rod.
     procedure :: downcast_assign => rod_downcast_assign

     !> Equality operator implementation.
     procedure :: rod_equals
     generic :: operator(==) => rod_equals

     !> Reads coordinates of a rod from string.
     procedure :: rod_from_str
     !> Writes a rod to a given Fortran output unit.
     procedure :: to_unit => rod_to_unit
     !> Serializes the coordinates of the rod to a JSON value.
     procedure :: coordinates_to_json => rod_coordinates_to_json
     !> Deserializes coordinates of the rod from a JSON value.
     procedure :: from_json => rod_from_json
     !> Returns the orientation vector of this rod.
     procedure :: orientation => rod_orientation
     !> Returns the orientation vector of this rod.
     procedure :: set_orientation => rod_set_orientation
     !> Translates the rod into a new position and rotates it to a new
     !! orientation.
     procedure :: move => rod_move
  end type rod

contains

  !> Reads @p this rod from @p str. If ios /= 0, the read failed.
  subroutine rod_from_str(this, str, ios)
    class(rod), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z, this%ux, &
         this%uy, this%uz
  end subroutine rod_from_str


  !> Returns the type of the rod in @p str.
  subroutine rod_typestr(str)
    character(len=:), allocatable, intent(out) :: str
    str = "rod"
  end subroutine rod_typestr

  
  !> Returns the descriptions of the rod coordinates in @p descr.
  subroutine rod_description(descr)
    character(kind=CK, len=3), allocatable, intent(inout) :: descr(:)
    descr = ["x  ", "y  ", "z  ", "ux ", "uy ", "uz "]
  end subroutine rod_description


  !> Downcasts @p a_particle to rod and assigns it to @p this.
  !! @p err == 3 if an error occurs.
  pure subroutine rod_downcast_assign(this, a_particle, err)
    class(rod), intent(inout) :: this
    class(particle), intent(in) :: a_particle
    integer, intent(out), optional :: err
    select type (a_particle)
    type is (rod)
       this = a_particle
       if (present(err)) err = 0
    class default
       if (present(err)) err = 3
    end select
  end subroutine rod_downcast_assign


  !> Assigns the rod @p another to @p this.
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


  !> Returns true if coordinates of @p another are equal to @p this.
  elemental function rod_equals(this, another) result(res)
    class(rod), intent(in) :: this
    type(rod), intent(in) :: another
    logical :: res
    res = (this%x == another%x) .and. (this%y == another%y) .and. &
         (this%z == another%z) .and. (this%ux == another%ux) .and. &
         (this%uy == another%uy) .and. (this%uz == another%uz)
  end function rod_equals


  !> Writes @p this rod to the Fortran output @p unit.
  subroutine rod_to_unit(this, unit)
    class(rod), intent(in) :: this
    integer, intent(in) :: unit
    write(unit, '(A,6(' // fmt_char_dp() // ',1X))', advance='no') &
         'rod', this%x, this%y, this%z, this%ux, this%uy, this%uz 
  end subroutine rod_to_unit


  !> Translates and rotates @p this to a new positiona and orientation.
  !! @p genstate is the random number generator to be used.
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


  !> Serializes the coordinates of @p this rod to the JSON value
  !! @p json_val.
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


  !> Deserializes the coordinates of @p this from the JSON value
  !! @p json_val.
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


  !> Returns the orientation vector of @p this.
  pure function rod_orientation(this) 
    class(rod), intent(in) :: this
    real(dp), dimension(3) :: rod_orientation
    rod_orientation = [this%ux, this%uy, this%uz]
  end function rod_orientation

  !> Sets the orientation of @p this rod to @p vec.
  pure subroutine rod_set_orientation(this, vec)
    class(rod), intent(inout) :: this
    real(dp), intent(in) :: vec(3)
    this%ux = vec(1)
    this%uy = vec(2)
    this%uz = vec(3)
  end subroutine rod_set_orientation
  
end module m_rod
