module m_particledat
  use iso_fortran_env, only: output_unit
  use m_particle, only: particle
  use num_kind, only: dp
  use utils, only: fmt_char_dp
  use particle_mover, only: transmove, rotate
  use json_module, only: json_value, CK, json_get, json_add, json_create_array
  use m_json_wrapper, only: get_parameter
  include 'rng.inc'
  implicit none
  
  !> Holds data of a uniaxial particle, e.g. Gay-Berne 
  !! particle. Using the rod variable one can use the same datatype to describe
  !! two different particles.
  !! 
  !! x, y, z are the cartesian positional coordinates of the particle.
  !! ux, uy, uz are the components of the unit vector that describes
  !! the 
  !! orientation of the particle.
  !! rod tells if the particle is a rodlike particle or e.g. spherically
  !! symmetric.
  !!
  !! :TODO: Remove rod variable. This module should only describe GB particles.
  !!
  type, extends(particle) :: particledat
     real(dp) :: ux = 0._dp
     real(dp) :: uy = 0._dp
     real(dp) :: uz = 1._dp
     logical :: rod = .true.
   contains
     procedure :: downcast_assign => particledat_downcast_assign
     procedure :: particledat_assign
     generic :: assignment(=) => particledat_assign
     procedure :: to_stdout
     procedure :: move => particledat_move
     procedure :: orientation => orientation
     procedure, nopass :: typestr => particledat_typestr
     procedure, nopass :: description => particledat_description
     procedure :: coordinates_to_json => particledat_coordinates_to_json
     procedure :: from_json => particledat_from_json
     procedure :: particledat_equals
     generic :: operator(==) => particledat_equals
  end type particledat
  
contains

  subroutine particledat_from_json(this, json_val)
    class(particledat), intent(inout) :: this
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
    call get_parameter(json_val, '[7]', this%rod)
  end subroutine particledat_from_json

  subroutine particledat_typestr(str)
    character(len=:), allocatable, intent(out) :: str
    str = "particledat"
  end subroutine particledat_typestr

  subroutine particledat_description(descr)
    character(kind=CK, len=3), allocatable, intent(inout) :: descr(:)
    descr = ["x  ", "y  ", "z  ", "ux ", "uy ", "uz ", "rod"]
  end subroutine particledat_description

  pure subroutine particledat_downcast_assign(this, a_particle, err)
    class(particledat), intent(inout) :: this
    class(particle), intent(in) :: a_particle
    integer, intent(out), optional :: err
    select type (a_particle)
    type is (particledat)
       this = a_particle
    class default
       if (present(err)) err = 3
    end select
  end subroutine particledat_downcast_assign

  elemental function particledat_equals(this, another) result(res)
    class(particledat), intent(in) :: this
    type(particledat), intent(in) :: another
    logical :: res
    res = (this%x == another%x) .and. (this%y == another%y) .and. &
         (this%z == another%z) .and. (this%ux == another%ux) .and. &
         (this%uy == another%uy) .and. (this%uz == another%uz) .and. &
         (this%rod .eqv. another%rod)
  end function particledat_equals
  
  !> Writes @p aparticle to the outputunit @p writeunit.
  subroutine writeparticle(writeunit, aparticle)
    integer, intent(in) :: writeunit
    type(particledat), intent(in) :: aparticle
    character(len = 2) :: typeid
    if (aparticle%rod) then
       typeid = 'gb'
    else 
       typeid = 'lj'
    end if
    write(writeunit, '(A,1X,6(' // fmt_char_dp() // ',1X))', advance='no') &
         typeid, aparticle%x, aparticle%y, aparticle%z, aparticle%ux, &
         aparticle%uy, aparticle%uz 
  end subroutine writeparticle
  
  !> Writes @p aparticle to the outputunit @p writeunit.
  subroutine to_stdout(this)
    class(particledat), intent(in) :: this
    character(len = 2) :: typeid
    if (this%rod) then
       typeid = 'gb'
    else 
       typeid = 'lj'
    end if
    write(output_unit, '(A,1X,6(' // fmt_char_dp() // ',1X))', advance='no') &
         typeid, this%x, this%y, this%z, this%ux, &
         this%uy, this%uz 
  end subroutine to_stdout

  subroutine particledat_coordinates_to_json(this, json_val)
  class(particledat), intent(in) :: this
  type(json_value), pointer :: json_val
  call json_create_array(json_val, '')
  call json_add(json_val, '', this%x)
  call json_add(json_val, '', this%y)
  call json_add(json_val, '', this%z)
  call json_add(json_val, '', this%ux)
  call json_add(json_val, '', this%uy)
  call json_add(json_val, '', this%uz)
  call json_add(json_val, '', this%rod)
end subroutine particledat_coordinates_to_json


!> Reads @p aparticle from @p readunit. If the reading succeeds, 
!! @p ios==0.
subroutine readparticle(readunit, aparticle, ios)
  type(particledat), intent(out) :: aparticle
  integer, intent(in) :: readunit
  integer, intent(out) :: ios
  character(len = 2) :: temp
  character(len=1000) :: buffer
  read(readunit, fmt='(A1000)', iostat=ios) buffer
  if (ios /= 0) return
  read(buffer, *, iostat=ios) temp
  if (temp == 'gb') then
     aparticle%rod = .true.
     read(buffer, fmt=*, iostat=ios) &
       temp, aparticle%x, aparticle%y, aparticle%z, aparticle%ux, &
       aparticle%uy, aparticle%uz
  else if (temp == 'lj') then 
     aparticle%rod = .false.
     read(buffer, fmt=*, iostat=ios) &
          temp, aparticle%x, aparticle%y, aparticle%z
  else
     ios = 999
  end if
  if (ios /= 0) backspace readunit
end subroutine readparticle

elemental subroutine particledat_assign(dest, src)
  class(particledat), intent(inout) :: dest
  type(particledat), intent(in) :: src
  dest%x = src%x
  dest%y = src%y
  dest%z = src%z
  dest%ux = src%ux
  dest%uy = src%uy
  dest%uz = src%uz
  dest%rod = src%rod
end subroutine particledat_assign

!> Returns the orientation of @p aparticle as a vector. @p vec must
!! be a unit vector.
pure function orientation(aparticle) 
  real(dp), dimension(3) :: orientation
  class(particledat), intent(in) :: aparticle
  orientation = (/aparticle%ux, aparticle%uy, aparticle%uz/)
end function orientation

!> Assigns the orientation @p vec to @p aparticle.
pure subroutine setorientation(aparticle, vec)
  type(particledat), intent(inout) :: aparticle
  real(dp), dimension(3), intent(in) :: vec
  aparticle%ux = vec(1)
  aparticle%uy = vec(2)
  aparticle%uz = vec(3)
end subroutine setorientation

!> Performs a combined translation and rotation of @p particle. 
!! @p genstate is the random number generator state.
pure subroutine particledat_move(this, genstate)    
  class(particledat), intent(inout) :: this
  type(rngstate), intent(inout) :: genstate
  real(dp) :: xn, yn, zn, uxn, uyn, uzn
  call transmove(this%x, this%y, this%z, xn, yn, zn, genstate)
  this%x = xn
  this%y = yn
  this%z = zn
  if (this%rod) then
     call rotate(this%ux, this%uy, this%uz, uxn, uyn, uzn, genstate)
     this%ux = uxn
     this%uy = uyn
     this%uz = uzn
  end if
end subroutine particledat_move


end module m_particledat
