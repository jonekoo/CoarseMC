module class_poly_box
use nrtype
use utils
implicit none

public :: poly_box
public :: new_box
public :: new_cylinder
public :: makebox
public :: createbox
public :: minimage
public :: mindistance
public :: gettype
public :: getx, gety, getz
public :: setx, sety, setz
public :: isxperiodic, isyperiodic, iszperiodic
public :: setxperiodicity, setyperiodicity, setzperiodicity
public :: volume
public :: writetostdio
public :: write
public :: binwrite
public :: read

private

character(len=15), dimension(2), parameter :: typeids = (/'rectangular', 'cylindrical'/)

!> Type to hold information about the geometry of the system.
!! Possible values of typeid are listed in the static array typeids.
!!
type poly_box
  private
  character(len=15) :: typeid = 'rectangular'
  real(dp) :: lx = 0._dp
  real(dp) :: ly = 0._dp
  real(dp) :: lz = 0._dp
  logical :: xperiodic = .true.
  logical :: yperiodic = .true.
  logical :: zperiodic = .true.
end type poly_box

interface new_box
  module procedure new_box1, new_box3
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

interface writetostdio
  module procedure pbox_writetostdio
end interface

interface volume
  module procedure pbox_volume
end interface

interface makebox
  module procedure pbox_makebox, pbox_copy
end interface

interface createbox
  module procedure pbox_createbox
end interface

interface write
  module procedure pbox_write
end interface

interface binwrite
  module procedure pbox_binwrite
end interface

interface read
  module procedure pbox_read
end interface

interface gettype
  module procedure pbox_gettype
end interface

contains 

  function new_box1(side) result(b)
    type(poly_box) :: b
    real(dp), intent(in) :: side
    b%lx = side
    b%ly = side
    b%lz = side
  end function

  function new_box3(lx, ly, lz) result(b)
    type(poly_box) :: b
    real(dp), intent(in) :: lx, ly, lz
    b%lx = lx
    b%ly = ly
    b%lz = lz
  end function

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
  end function

  !! Returns a "well-defined" periodic cubic box with sidelengths @p side.
  !! This routine is meant to be a initializing routine after which 
  !! customizations can be done.
  !!
  !! @p side the length of one side.
  !!
  !! @return a cubic box.
  !!  
  subroutine pbox_makebox(simbox, side) 
    type(poly_box), intent(out) :: simbox
    real(dp), intent(in) :: side
    simbox%lx = side
    simbox%ly = side
    simbox%lz = side
  end subroutine

  subroutine pbox_copy(copy, simbox)
    type(poly_box), intent(out) :: copy
    type(poly_box), intent(in) :: simbox
    copy = simbox
  end subroutine

  subroutine pbox_createbox(bp, boxstring)
    type(poly_box), intent(inout) :: bp
    character(len = *), intent(in) :: boxstring
    character(len = 5) :: typeid
    integer :: ios
    read(boxstring, *, iostat = ios) typeid, bp%lx, bp%ly, bp%lz, bp%xperiodic, bp%yperiodic, bp%zperiodic
    if (ios /= 0) then
      write(6, *) 'pbox_createbox: failed to read box with iostat = ', ios, '. Stopping.'
      stop 
    else if('box' /= typeid) then
      write(6, *) 'pbox_createbox: Warning converting from ', trim(adjustl(typeid)), 'to box!'
    end if
  end subroutine

  subroutine pbox_write(unit, bp)
    type(poly_box), intent(in) :: bp
    integer, intent(in) :: unit
    write(unit, '(A15, 3' // fmt_char_dp() // ',3L3)', advance='no') bp%typeid, &
    bp%lx, bp%ly, bp%lz, bp%xperiodic, bp%yperiodic, bp%zperiodic
    !write(line, '(A, 3' // fmt_char_dp() // ')') typeid, bp%lx, bp%ly, bp%lz
    !write(*, *) trim(line)
  end subroutine

  subroutine pbox_binwrite(unit, bp)
    type(poly_box), intent(in) :: bp
    integer, intent(in) :: unit
    write(unit) bp%typeid, bp%lx, bp%ly, bp%lz, bp%xperiodic, bp%yperiodic,&
    bp%zperiodic
  end subroutine

  function pbox_istypeid(string) result(isvalid)
    character(len=*), intent(in) :: string
    logical :: isvalid
    if (any(typeids == string)) then
      isvalid = .true.
    else
      isvalid = .false.
    end if
  end function

  subroutine pbox_read(unit, bp, ios)
    integer, intent(in) :: unit
    type(poly_box), intent(out) :: bp
    integer :: ios
    read(unit, *, iostat=ios) bp%typeid, bp%lx, bp%ly, bp%lz, bp%xperiodic, &
    bp%yperiodic, bp%zperiodic
    if (.not. pbox_istypeid(bp%typeid)) ios = 999
    if (ios /= 0) then
      !! Return unit to its previous state if reading fails.
      backspace unit
    end if
  end subroutine

  elemental function pbox_gettype(abox)
    type(poly_box), intent(in) :: abox
    character(len=len(abox%typeid)) pbox_gettype
    pbox_gettype = abox%typeid
  end function

  real(dp) elemental function pbox_getx(abox)
    type(poly_box), intent(in) :: abox
    pbox_getx = abox%lx
  end function

  real(dp) elemental function pbox_gety(abox)
    type(poly_box), intent(in) :: abox
    pbox_gety = abox%ly
  end function

  real(dp) elemental function pbox_getz(abox)
    type(poly_box), intent(in) :: abox
    pbox_getz = abox%lz
  end function

  elemental function pbox_isxperiodic(abox) result(isperiodic)
    type(poly_box), intent(in) :: abox
    logical :: isperiodic
    isperiodic = abox%xperiodic
  end function

  elemental function pbox_isyperiodic(abox) result(isperiodic)
    type(poly_box), intent(in) :: abox
    logical :: isperiodic
    isperiodic = abox%yperiodic
  end function

  elemental function pbox_iszperiodic(abox) result(isperiodic)
    type(poly_box), intent(in) :: abox
    logical :: isperiodic
    isperiodic = abox%zperiodic
  end function

  pure subroutine pbox_setx(abox, x)
    implicit none
    type(poly_box), intent(inout) :: abox
    real(dp), intent(in) :: x
    abox%lx = x
  end subroutine 
  
  pure subroutine pbox_sety(abox, y)
    implicit none
    type(poly_box), intent(inout) :: abox
    real(dp), intent(in) :: y
    abox%ly = y
  end subroutine 
  
  pure subroutine pbox_setz(abox, z)
    implicit none
    type(poly_box), intent(inout) :: abox
    real(dp), intent(in) :: z
    abox%lz = z
  end subroutine 

  pure subroutine pbox_setxperiodicity(abox, isperiodic)
    type(poly_box), intent(inout) :: abox
    logical, intent(in) :: isperiodic
    abox%xperiodic = isperiodic
  end subroutine

  pure subroutine pbox_setyperiodicity(abox, isperiodic)
    type(poly_box), intent(inout) :: abox
    logical, intent(in) :: isperiodic
    abox%yperiodic = isperiodic
  end subroutine

  pure subroutine pbox_setzperiodicity(abox, isperiodic)
    type(poly_box), intent(inout) :: abox
    logical, intent(in) :: isperiodic
    abox%zperiodic = isperiodic
  end subroutine

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
  end function

  pure function pbox_mindistance(simbox, ri, rj)
    real(dp) :: pbox_mindistance
    type(poly_box), intent(in) :: simbox
    real(dp), dimension(3), intent(in) :: ri
    real(dp), dimension(3), intent(in) :: rj
    real(dp), dimension(3) :: rij
    rij = minimage(simbox, rj - ri)
    pbox_mindistance = sqrt(dot_product(rij, rij))
  end function 

  pure recursive function pbox_minimage(simbox, r)
    implicit none
    real(dp), dimension(3) :: pbox_minimage
    type(poly_box), intent(in) :: simbox
    real(dp), dimension(3), intent(in) :: r
    !! Make periodic transformations
    pbox_minimage = r
    if (simbox%xperiodic) then
      !pbox_minimage(1) = r(1) - simbox%lx * anint(r(1)/simbox%lx)
      if (pbox_minimage(1) < -0.5_dp*simbox%lx) then
        pbox_minimage(1)=pbox_minimage(1)+simbox%lx 
      else if(pbox_minimage(1) >= 0.5_dp*simbox%lx) then
        pbox_minimage(1)=pbox_minimage(1)-simbox%lx
      end if
    end if
    if (simbox%yperiodic) then
      !pbox_minimage(2) = r(2) - simbox%ly * anint(r(2)/simbox%ly)
      if (pbox_minimage(2) < -0.5_dp*simbox%ly) then
        pbox_minimage(2)=pbox_minimage(2)+simbox%ly 
      else if(pbox_minimage(2) >= 0.5_dp*simbox%ly) then
        pbox_minimage(2)=pbox_minimage(2)-simbox%ly
      end if
    end if
    if (simbox%zperiodic) then
      !pbox_minimage(3) = r(3) - simbox%lz * anint(r(3)/simbox%lz)
      !! This is faster than the solution above:
      if (pbox_minimage(3) < -0.5_dp*simbox%lz) then
        pbox_minimage(3)=pbox_minimage(3)+simbox%lz 
      else if(pbox_minimage(3) >= 0.5_dp*simbox%lz) then
        pbox_minimage(3)=pbox_minimage(3)-simbox%lz
      end if
    end if
    ! This seems to be slower when optimized:
    !pbox_minimage(1) = r(1) - simbox%lx * count((/simbox%xperiodic/)) * anint(r(1)/simbox%lx) 
    !pbox_minimage(2) = r(2) - simbox%ly * count((/simbox%yperiodic/)) * anint(r(2)/simbox%ly) 
    !pbox_minimage(3) = r(3) - simbox%lz * count((/simbox%zperiodic/)) * anint(r(3)/simbox%lz) 
  end function

  subroutine pbox_writetostdio(simbox)
    type(poly_box), intent(in) :: simbox
    write(*,*) simbox
  end subroutine

end module
