!> Contains the doublerod type and related procedures.
!> forked from m_rod by Perttu Tuovinen
module m_doublerod
  use iso_fortran_env, only: output_unit
  use m_point, only: point
  use num_kind, only: dp
  use utils, only: fmt_char_dp
  use particle_mover, only: transmove, rotate, doublerotate
  use json_module, only: json_create_array, CK, json_value, json_add
  use m_json_wrapper, only: get_parameter
  use mt_stream, only: rng=>genrand_double1_s, rngstate=>mt_state
  implicit none

  !> A particle consisting of two rods attached to each other on one end.
  type, extends(point) :: doublerod
     !> The components of the unit vectors of orientation.
     real(dp) :: ux = 0._dp
     real(dp) :: uy = 0._dp
     real(dp) :: uz = 1._dp
     real(dp) :: vx = 0._dp
     real(dp) :: vy = 0._dp
     real(dp) :: vz = -1._dp     
   contains
     !> Returns the type of the double rod as a string.
     procedure, nopass :: typestr => doublerod_typestr
     !> Returns the description of the coordinates for this doublerod as an
     !! array of strings.
     procedure, nopass :: description => doublerod_description

     !> Assignment operator implementation.
     procedure :: doublerod_assign
     generic :: assignment(=) => doublerod_assign
     !> Downcasts and assigns a particle to a doublerod.
     procedure :: downcast_assign => doublerod_downcast_assign

     !> Equality operator implementation.
     procedure :: doublerod_equals
     generic :: operator(==) => doublerod_equals

     !> Reads coordinates of a doublerod from string.
     procedure :: doublerod_from_str
     !> Writes a doublerod to a given Fortran output unit.
     procedure :: to_unit => doublerod_to_unit
     !> Serializes the coordinates of the doublerod to a JSON value.
     procedure :: coordinates_to_json => doublerod_coordinates_to_json
     !> Deserializes coordinates of the doublerod from a JSON value.
     procedure :: from_json => doublerod_from_json



     !> Return the orientation vectors of this doublerod.
     procedure :: orientation => doublerod_orientation_u
     procedure :: orientation_v => doublerod_orientation_v
     procedure :: orientation_u => doublerod_orientation_u

     !> Set the orientation vectors of this doublerod.
     procedure :: set_orientation => doublerod_set_orientation
     procedure :: set_orientation_u => doublerod_set_orientation_u
     procedure :: set_orientation_v => doublerod_set_orientation_v


     !> Translates the doublerod into a new position and rotates it to a new
     !! orientation.
     procedure :: move => doublerod_move
  end type doublerod

contains

  !> Reads @p this doublerod from @p str. If ios /= 0, the read failed.
  subroutine doublerod_from_str(this, str, ios)
    class(doublerod), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, intent(out) :: ios
    read(str, fmt=*, iostat=ios) this%x, this%y, this%z, this%ux, &
         this%uy, this%uz, this%vx, this%vy, this%vz
  end subroutine doublerod_from_str


  !> Returns the type of the doublerod in @p str.
  subroutine doublerod_typestr(str)
    character(len=:), allocatable, intent(out) :: str
    str = "doublerod"
  end subroutine doublerod_typestr

  
  !> Returns the descriptions of the doublerod coordinates in @p descr.
  subroutine doublerod_description(descr)
    character(kind=CK, len=3), allocatable, intent(inout) :: descr(:)
    descr = ["x  ", "y  ", "z  ", "ux ", "uy ", "uz ", "vx ", "vy ", "vz "]
  end subroutine doublerod_description


  !> Downcasts @p a_particle to doublerod and assigns it to @p this.
  !! @p err == 3 if an error occurs.
  pure subroutine doublerod_downcast_assign(this, a_particle, err)
    class(doublerod), intent(inout) :: this
    class(point), intent(in) :: a_particle
    integer, intent(out), optional :: err
    select type (a_particle)
    type is (doublerod)
       this = a_particle
       if (present(err)) err = 0
    class default
       if (present(err)) err = 6
       ! this "6" is just made up on the fly lol
    end select
  end subroutine doublerod_downcast_assign


  !> Assigns the doublerod @p another to @p this.
  pure subroutine doublerod_assign(this, another)
    class(doublerod), intent(inout) :: this
    type(doublerod), intent(in) :: another
    this%x = another%x
    this%y = another%y
    this%z = another%z
    this%ux = another%ux
    this%uy = another%uy
    this%uz = another%uz
    this%vx = another%vx
    this%vy = another%vy
    this%vz = another%vz
  end subroutine doublerod_assign


  !> Returns true if coordinates of @p another are equal to @p this.
  elemental function doublerod_equals(this, another) result(res)
    class(doublerod), intent(in) :: this
    type(doublerod), intent(in) :: another
    logical :: res
    res = (this%x == another%x) .and. (this%y == another%y) .and. &
         (this%z == another%z) .and. (this%ux == another%ux) .and. &
         (this%uy == another%uy) .and. (this%uz == another%uz) .and. &
         (this%vx == another%vx) .and. (this%vy == another%vy) .and. &
         (this%vz == another%vz)
  end function doublerod_equals


  !> Writes @p this doublerod to the Fortran output @p unit.
  subroutine doublerod_to_unit(this, unit)
    class(doublerod), intent(in) :: this
    integer, intent(in) :: unit
    write(unit, '(A,6(' // fmt_char_dp() // ',1X))', advance='no') &
         'doublerod', this%x, this%y, this%z, this%ux, this%uy, this%uz, &
         this%vx, this%vy, this%vz 
  end subroutine doublerod_to_unit


  !> Translates and rotates @p this to a new position and orientation.
  !! @p genstate is the random number generator to be used.
  ! CHECK if this -really- rotates vector v correctly
  pure subroutine doublerod_move(this, genstate)    
    class(doublerod), intent(inout) :: this
    type(rngstate), intent(inout) :: genstate
    real(dp) :: xn, yn, zn, uxn, uyn, uzn, vxn, vyn, vzn
    call transmove(this%x, this%y, this%z, xn, yn, zn, genstate)
    this%x = xn
    this%y = yn
    this%z = zn
    call doublerotate(this%ux, this%uy, this%uz, this%vx, this%vy, &
      this%vz, uxn, uyn, uzn, vxn, vyn, vzn, genstate)
    this%ux = uxn 
    this%uy = uyn
    this%uz = uzn
    this%vx = vxn
    this%vy = vyn
    this%vz = vzn
  end subroutine doublerod_move


  !> Serializes the coordinates of @p this doublerod to the JSON value
  !! @p json_val.
  subroutine doublerod_coordinates_to_json(this, json_val)
    class(doublerod), intent(in) :: this
    type(json_value), pointer :: json_val
    call json_create_array(json_val, '')
    call json_add(json_val, '', this%x)
    call json_add(json_val, '', this%y)
    call json_add(json_val, '', this%z)
    call json_add(json_val, '', this%ux)
    call json_add(json_val, '', this%uy)
    call json_add(json_val, '', this%uz)
    call json_add(json_val, '', this%vx)
    call json_add(json_val, '', this%vy)
    call json_add(json_val, '', this%vz)
  end subroutine doublerod_coordinates_to_json


  !> Deserializes the coordinates of @p this from the JSON value
  !! @p json_val.
  subroutine doublerod_from_json(this, json_val)
    class(doublerod), intent(inout) :: this
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
    call get_parameter(json_val, '[7]', this%vx, error_lb=-1._dp, &
         error_ub=1._dp)
    call get_parameter(json_val, '[8]', this%vy, error_lb=-1._dp, &
         error_ub=1._dp)
    call get_parameter(json_val, '[9]', this%vz, error_lb=-1._dp, &
         error_ub=1._dp)
  end subroutine doublerod_from_json


!yep these be needed

  !> Returns the orientation vector u of @p this.
  ! this and the following v counterpart aren't inherited but exclusive to doublerod.
  pure function doublerod_orientation_u(this) 
    class(doublerod), intent(in) :: this
    real(dp), dimension(3) :: doublerod_orientation_u
    doublerod_orientation_u = [this%ux, this%uy, this%uz]
  end function doublerod_orientation_u

  !> Returns the orientation vector v of @p this.
  pure function doublerod_orientation_v(this) 
    class(doublerod), intent(in) :: this
    real(dp), dimension(3) :: doublerod_orientation_v
    doublerod_orientation_v = [this%vx, this%vy, this%vz]
  end function doublerod_orientation_v

  !> Sets the orientation of @p this doublerod u and v vectors to @p vec 
  !! and its inverse, respectively. Thus, this works exclusively for 
  !! trans-azobenzene type molecules.
  pure subroutine doublerod_set_orientation(this, vec)
    class(doublerod), intent(inout) :: this
    real(dp), intent(in) :: vec(3)
    this%ux = vec(1)
    this%uy = vec(2)
    this%uz = vec(3)
    this%vx = -vec(1)
    this%vy = -vec(2)
    this%vz = -vec(3)
  end subroutine doublerod_set_orientation

  !> Sets the orientation of @p this doublerod u vector to @p vec.
  pure subroutine doublerod_set_orientation_u(this, vec)
    class(doublerod), intent(inout) :: this
    real(dp), intent(in) :: vec(3)
    this%ux = vec(1)
    this%uy = vec(2)
    this%uz = vec(3)
  end subroutine doublerod_set_orientation_u

   !> Sets the orientation of @p this doublerod v vector to @p vec.
  pure subroutine doublerod_set_orientation_v(this, vec)
    class(doublerod), intent(inout) :: this
    real(dp), intent(in) :: vec(3)
    this%vx = vec(1)
    this%vy = vec(2)
    this%vz = vec(3)
  end subroutine doublerod_set_orientation_v
  
end module m_doublerod
