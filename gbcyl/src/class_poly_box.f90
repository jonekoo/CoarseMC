!> Implements routines related to the simulation cell/box, such as
!! minimum image calculations.
module class_poly_box
use nrtype
use utils
use iso_fortran_env
implicit none

public :: poly_box
public :: new_box
public :: new_cylinder
public :: minimage
public :: mindistance
public :: gettype
public :: getx, gety, getz
public :: setx, sety, setz
public :: isxperiodic, isyperiodic, iszperiodic
public :: setxperiodicity, setyperiodicity, setzperiodicity
public :: volume
public :: pbox_read
public :: pbox_write

private

character(len=15), dimension(2), parameter :: typeids = (/'rectangular', 'cylindrical'/)

!> Stores the geometry of the system. Possible values of typeid are
!! listed in the static array typeids.
type poly_box
  character(len=15) :: typeid = 'rectangular'
  real(dp) :: lx = 0._dp
  real(dp) :: ly = 0._dp
  real(dp) :: lz = 0._dp
  logical :: xperiodic = .true.
  logical :: yperiodic = .true.
  logical :: zperiodic = .true.
end type poly_box

interface new_box
  module procedure new_box_all, new_box1, new_box3
end interface

interface minimage
  module procedure pbox_minimage
end interface

interface mindistance
  module procedure pbox_mindistance
end interface

interface setx
  module procedure pbox_setx
end interface

interface sety
  module procedure pbox_sety
end interface

interface setz
  module procedure pbox_setz
end interface

interface getx
  module procedure pbox_getx
end interface

interface gety
  module procedure pbox_gety
end interface

interface getz
  module procedure pbox_getz
end interface

interface isxperiodic
  module procedure pbox_isxperiodic
end interface

interface isyperiodic
  module procedure pbox_isyperiodic
end interface

interface iszperiodic
  module procedure pbox_iszperiodic
end interface

interface setxperiodicity
  module procedure pbox_setxperiodicity
end interface

interface setyperiodicity
  module procedure pbox_setyperiodicity
end interface

interface setzperiodicity
  module procedure pbox_setzperiodicity
end interface

interface volume
  module procedure pbox_volume
end interface

interface gettype
  module procedure pbox_gettype
end interface

contains 

!> Returns a poly_box with the given typeid @p boxtype. @p lx and, 
!! optionally, @p ly and @p lz give the dimensions of the box. If
!! @p boxtype == 'rectangular' all three dimensions are considered. If
!! @p boxtype == 'cylinder', the diameter of the cylinder is set as the
!! greater of @p lx and @p ly. @p lz is then the height of the cylinder.
!! When @p ly or @p lz is not given, it is set equal to @p lx.
function new_box_all(boxtype, lx, ly, lz) result(b)
  character(len=*), intent(in) :: boxtype
  real(dp), intent(in) :: lx
  real(dp), intent(in), optional :: ly, lz
  real(dp) :: ly_, lz_ 
  type(poly_box) :: b
  ly_ = lx
  lz_ = lx
  if (present(ly)) ly_ = ly
  if (present(lz)) lz_ = lz
  if (boxtype == 'rectangular') then
     b = new_box3(lx, ly_, lz_)        
  else if(boxtype == 'cylindrical') then
     b = new_cylinder(max(lx, ly_), lz_)
  else if (.not. pbox_istypeid(boxtype)) then
     write(error_unit, *) 'Error: ', boxtype, 'is not a valid typeid for' // &
          ' poly_box.'
     write(error_unit, *) 'Valid options are', typeids(:)
     stop
  else
     write(error_unit, *) 'Error: class_poly_box:new_box_all: internal ' // &
          'consistency error. Please send a bug report to the developer!'
  end if
end function

!> Returns a rectangular type poly_box with sides of equal length
!! @p side and periodic boundaries in all directions.
function new_box1(side) result(b)
  type(poly_box) :: b
  real(dp), intent(in) :: side
  b%lx = side
  b%ly = side
  b%lz = side
end function new_box1

!> Returns a rectangular type poly_box with sides of length @p lx,ly,lz
!! and periodic boundaries in all directions.
function new_box3(lx, ly, lz) result(b)
  type(poly_box) :: b
  real(dp), intent(in) :: lx, ly, lz
  b%lx = lx
  b%ly = ly
  b%lz = lz
end function new_box3

!> Returns a cylindrical poly_box with @p diameter and @p height. The
!! height of the cylinder is in z-direction. z-direction is periodic. 
function new_cylinder(diameter, height) result(pbox) 
  real(dp), intent(in) :: diameter
  real(dp), intent(in) :: height
  type(poly_box) :: pbox
  pbox%typeid = typeids(2)
  pbox%lx = diameter
  pbox%ly = diameter
  pbox%lz = height
  pbox%xperiodic = .false.
  pbox%yperiodic = .false.
  pbox%zperiodic = .true.
end function new_cylinder

!> Writes the poly_box @p bp to the given output @p unit.
subroutine pbox_write(unit, bp)
  type(poly_box), intent(in) :: bp
  integer, intent(in) :: unit
  write(unit, '(A15, 3' // fmt_char_dp() // ',3L3)', advance='no') bp%typeid, &
       bp%lx, bp%ly, bp%lz, bp%xperiodic, bp%yperiodic, bp%zperiodic
end subroutine pbox_write

!> Returns true if typeid @p string is valid typeid for a poly_box.
function pbox_istypeid(string) result(isvalid)
  character(len=*), intent(in) :: string
  logical :: isvalid
  if (any(typeids == string)) then
     isvalid = .true.
  else
     isvalid = .false.
  end if
end function pbox_istypeid

!> Reads the poly_box @p bp from @p unit. If the read fails, 
!! @p ios /= 0.
subroutine pbox_read(unit, bp, ios)
  integer, intent(in) :: unit
  type(poly_box), intent(inout) :: bp
  integer :: ios
  read(unit, *, iostat=ios) bp%typeid, bp%lx, bp%ly, bp%lz, bp%xperiodic, &
       bp%yperiodic, bp%zperiodic
  if (.not. pbox_istypeid(bp%typeid)) ios = 999
  if (ios /= 0) then
     !! Return unit to its previous state if reading fails.
     backspace unit
  end if
end subroutine pbox_read

!> Returns the typeid of the @p abox in a string.
elemental function pbox_gettype(abox)
  type(poly_box), intent(in) :: abox
  character(len=len(abox%typeid)) pbox_gettype
  pbox_gettype = abox%typeid
end function pbox_gettype

!> Returns the size of @p abox in the x-direction. 
real(dp) elemental function pbox_getx(abox)
  type(poly_box), intent(in) :: abox
  pbox_getx = abox%lx
end function pbox_getx

!> Returns the size of @p abox in the y-direction. 
real(dp) elemental function pbox_gety(abox)
  type(poly_box), intent(in) :: abox
  pbox_gety = abox%ly
end function pbox_gety

!> Returns the size of @p abox in the z-direction. 
real(dp) elemental function pbox_getz(abox)
  type(poly_box), intent(in) :: abox
  pbox_getz = abox%lz
end function pbox_getz

!> Returns true if @p abox is periodic in the x-direction.
elemental function pbox_isxperiodic(abox) result(isperiodic)
  type(poly_box), intent(in) :: abox
  logical :: isperiodic
  isperiodic = abox%xperiodic
end function pbox_isxperiodic

!> Returns true if @p abox is periodic in the y-direction.
elemental function pbox_isyperiodic(abox) result(isperiodic)
  type(poly_box), intent(in) :: abox
  logical :: isperiodic
  isperiodic = abox%yperiodic
end function pbox_isyperiodic

!> Returns true if @p abox is periodic in the z-direction.
elemental function pbox_iszperiodic(abox) result(isperiodic)
  type(poly_box), intent(in) :: abox
  logical :: isperiodic
  isperiodic = abox%zperiodic
end function pbox_iszperiodic

!> Sets the length of @p abox to @p x in the x-direction.
pure subroutine pbox_setx(abox, x)
  type(poly_box), intent(inout) :: abox
  real(dp), intent(in) :: x
  abox%lx = x
  if (abox%typeid == 'cylindrical') abox%ly = x
end subroutine pbox_setx
  
!> Sets the length of @p abox to @p y in the y-direction.
pure subroutine pbox_sety(abox, y)
  type(poly_box), intent(inout) :: abox
  real(dp), intent(in) :: y
  abox%ly = y
end subroutine pbox_sety

!> Sets the length of @p abox to @p z in the z-direction.
pure subroutine pbox_setz(abox, z)
  type(poly_box), intent(inout) :: abox
  real(dp), intent(in) :: z
  abox%lz = z
end subroutine pbox_setz

!> Sets the periodicity of @p abox to @p isperiodic in the
!! x-direction.
pure subroutine pbox_setxperiodicity(abox, isperiodic)
  type(poly_box), intent(inout) :: abox
  logical, intent(in) :: isperiodic
  abox%xperiodic = isperiodic
end subroutine pbox_setxperiodicity

!> Sets the periodicity of @p abox to @p isperiodic in the
!! y-direction.
pure subroutine pbox_setyperiodicity(abox, isperiodic)
  type(poly_box), intent(inout) :: abox
  logical, intent(in) :: isperiodic
  abox%yperiodic = isperiodic
end subroutine pbox_setyperiodicity

!> Sets the periodicity of @p abox to @p isperiodic in the
!! z-direction.
pure subroutine pbox_setzperiodicity(abox, isperiodic)
  type(poly_box), intent(inout) :: abox
  logical, intent(in) :: isperiodic
  abox%zperiodic = isperiodic
end subroutine pbox_setzperiodicity

!> Returns the volume of @p simbox.
pure function pbox_volume(simbox)
  real(dp) :: pbox_volume
  type(poly_box), intent(in) :: simbox
  if (simbox%typeid == 'rectangular') then
     pbox_volume = simbox%lx * simbox%ly * simbox%lz
  else if (simbox%typeid == 'cylindrical') then
     pbox_volume = atan(1._dp) * simbox%lx**2 * simbox%lz
     !! The above translates to pi/4*diameter^2*lz = pi*(diameter/2 ^2*lz
  else 
     pbox_volume = -1._dp
  end if
end function pbox_volume

!> Returns the minimum image distance between points @p ri and @p rj in
!! @p simbox.
pure function pbox_mindistance(simbox, ri, rj)
  real(dp) :: pbox_mindistance
  type(poly_box), intent(in) :: simbox
  real(dp), dimension(3), intent(in) :: ri
  real(dp), dimension(3), intent(in) :: rj
  real(dp), dimension(3) :: rij
  rij = minimage(simbox, rj - ri)
  pbox_mindistance = sqrt(dot_product(rij, rij))
end function pbox_mindistance

!> Returns the minimum image vector in @p simbox corresponding to @p r.
pure function pbox_minimage(simbox, r) result(mi)
  real(dp), dimension(3) :: mi
  type(poly_box), intent(in) :: simbox
  real(dp), dimension(3), intent(in) :: r
  !! Make periodic transformations
  mi = r
  if (simbox%xperiodic) then
     if (mi(1) < -0.5 * simbox%lx) then
        mi(1) = mi(1)+simbox%lx 
     else if(mi(1) >= 0.5 * simbox%lx) then
        mi(1) = mi(1) - simbox%lx
     end if
  end if
  if (simbox%yperiodic) then
     if (mi(2) < -0.5 * simbox%ly) then
        mi(2) = mi(2) + simbox%ly 
     else if(mi(2) >= 0.5 * simbox%ly) then
        mi(2) = mi(2) - simbox%ly
     end if
  end if
  if (simbox%zperiodic) then
     !! This is faster than the solution above:
     if (mi(3) < -0.5 * simbox%lz) then
        mi(3) = mi(3) + simbox%lz 
     else if(mi(3) >= 0.5 * simbox%lz) then
        mi(3) = mi(3) - simbox%lz
     end if
  end if
end function pbox_minimage

end module
