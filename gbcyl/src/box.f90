module box
use nrtype
use particle
implicit none

public :: boxdat
public :: new_box
public :: makebox
public :: createbox
public :: minimage
public :: mindistance
public :: getx, gety, getz
public :: setx, sety, setz
public :: isxperiodic, isyperiodic, iszperiodic
public :: setxperiodicity, setyperiodicity, setzperiodicity
public :: volume
public :: scale
public :: writetostdio

private

type boxdat
  private
  real(dp) :: lx
  real(dp) :: ly
  real(dp) :: lz
  logical :: xperiodic
  logical :: yperiodic
  logical :: zperiodic    
end type boxdat

interface minimage
  module procedure minimage1, minimage2
end interface

interface mindistance
  module procedure box_mindistance
end interface

interface scale
  module procedure box_scale
end interface

interface setx
  module procedure box_setx
end interface

interface sety
  module procedure box_sety
end interface

interface setz
  module procedure box_setz
end interface

interface getx
  module procedure box_getx
end interface

interface gety
  module procedure box_gety
end interface

interface getz
  module procedure box_getz
end interface

interface isxperiodic
  module procedure box_isxp
end interface

interface isyperiodic
  module procedure box_isyperiodic
end interface

interface iszperiodic
  module procedure box_iszperiodic
end interface

interface setxperiodicity
  module procedure box_setxperiodicity
end interface

interface setyperiodicity
  module procedure box_setyperiodicity
end interface

interface setzperiodicity
  module procedure box_setzperiodicity
end interface

interface writetostdio
  module procedure box_writetostdio
end interface

interface volume
  module procedure box_volume
end interface

interface makebox
  module procedure box_makebox, box_copy
end interface

interface createbox
  module procedure box_createbox
end interface

contains 

  !! :TODO: Assignment operator, if needed. 

  function new_box(side) result(b)
    type(boxdat) :: b
    real(dp), intent(in) :: side
    call box_makebox(b, side)
  end function

  !! Returns a "well-defined" periodic cubic box with sidelengths @p side.
  !! This routine is meant to be a initializing routine after which 
  !! customizations can be done.
  !!
  !! @p side the length of one side.
  !!
  !! @return a cubic box.
  !!  
  subroutine box_makebox(simbox, side) 
    type(boxdat), intent(out) :: simbox
    real(dp), intent(in) :: side
    simbox%lx = side
    simbox%ly = side
    simbox%lz = side
    simbox%xperiodic = .true.
    simbox%yperiodic = .true.
    simbox%zperiodic = .true.
  end subroutine

  subroutine box_copy(copy, simbox)
    type(boxdat), intent(out) :: copy
    type(boxdat), intent(in) :: simbox
    copy = simbox
  end subroutine

  subroutine box_createbox(bp, boxstring)
    type(boxdat), intent(inout) :: bp
    character(len = *), intent(in) :: boxstring
    character(len = 5) :: typeid
    integer :: ios
    read(boxstring, *, iostat = ios) typeid, bp%lx, bp%ly, bp%lz, bp%xperiodic, bp%yperiodic, bp%zperiodic
    if (ios /= 0) then
      write(6, *) 'box_createbox: failed to read box with iostat = ', ios, '. Stopping.'
      stop 
    else if('box' /= typeid) then
      write(6, *) 'box_createbox: Warning converting from ', trim(adjustl(typeid)), 'to box!'
    end if
  end subroutine

  real(dp) elemental function box_getx(abox)
    type(boxdat), intent(in) :: abox
    box_getx = abox%lx
  end function

  real(dp) elemental function box_gety(abox)
    type(boxdat), intent(in) :: abox
    box_gety = abox%ly
  end function

  real(dp) elemental function box_getz(abox)
    type(boxdat), intent(in) :: abox
    box_getz = abox%lz
  end function

  elemental function box_isxp(abox) result(isperiodic)
    type(boxdat), intent(in) :: abox
    logical :: isperiodic
    isperiodic = abox%xperiodic
  end function

  elemental function box_isyperiodic(abox) result(isperiodic)
    type(boxdat), intent(in) :: abox
    logical :: isperiodic
    isperiodic = abox%yperiodic
  end function

  elemental function box_iszperiodic(abox) result(isperiodic)
    type(boxdat), intent(in) :: abox
    logical :: isperiodic
    isperiodic = abox%zperiodic
  end function

  pure subroutine box_setx(abox, x)
    implicit none
    type(boxdat), intent(inout) :: abox
    real(dp), intent(in) :: x
    abox%lx = x
  end subroutine 
  
  pure subroutine box_sety(abox, y)
    implicit none
    type(boxdat), intent(inout) :: abox
    real(dp), intent(in) :: y
    abox%ly = y
  end subroutine 
  
  pure subroutine box_setz(abox, z)
    implicit none
    type(boxdat), intent(inout) :: abox
    real(dp), intent(in) :: z
    abox%lz = z
  end subroutine 

  pure subroutine box_setxperiodicity(abox, isperiodic)
    type(boxdat), intent(inout) :: abox
    logical, intent(in) :: isperiodic
    abox%xperiodic = isperiodic
  end subroutine

  pure subroutine box_setyperiodicity(abox, isperiodic)
    type(boxdat), intent(inout) :: abox
    logical, intent(in) :: isperiodic
    abox%yperiodic = isperiodic
  end subroutine

  pure subroutine box_setzperiodicity(abox, isperiodic)
    type(boxdat), intent(inout) :: abox
    logical, intent(in) :: isperiodic
    abox%zperiodic = isperiodic
  end subroutine

  pure function box_volume(simbox)
    real(dp) :: box_volume
    type(boxdat), intent(in) :: simbox
    box_volume = simbox%lx * simbox%ly * simbox%lz
  end function

  pure function box_mindistance(simbox, ri, rj)
    real(dp) :: box_mindistance
    type(boxdat), intent(in) :: simbox
    real(dp), dimension(3), intent(in) :: ri
    real(dp), dimension(3), intent(in) :: rj
    real(dp), dimension(3) :: rij
    rij = minimage(simbox, ri, rj)
    box_mindistance = sqrt(dot_product(rij, rij))
  end function 

  pure function minimage2(simbox, ri, rj) result(rij)
    real(dp), dimension(3) :: rij
    type(boxdat), intent(in) :: simbox
    real(dp), dimension(3), intent(in) :: ri
    real(dp), dimension(3), intent(in) :: rj
    rij = rj - ri
    rij = minimage1(simbox, rij)
  end function

  pure recursive function minimage1(simbox, r)
    implicit none
    real(dp), dimension(3) :: minimage1
    type(boxdat), intent(in) :: simbox
    real(dp), dimension(3), intent(in) :: r
    !! Make periodic transformations
    minimage1 = r
    if (simbox%xperiodic) then
      minimage1(1) = r(1) - simbox%lx * anint(r(1)/simbox%lx)
    end if
    if (simbox%yperiodic) then
      minimage1(2) = r(2) - simbox%ly * anint(r(2)/simbox%ly)
    end if
    if (simbox%zperiodic) then
      minimage1(3) = r(3) - simbox%lz * anint(r(3)/simbox%lz)
    end if
  end function

  function box_scale(simbox, maxscaling, rng) result(scaling)
    type(boxdat), intent(inout) :: simbox
    real(dp), intent(in) :: maxscaling
    real(dp) :: dx, dy, dz
    interface
      function rng()
        real(8) :: rng
      end function
    end interface
    real(dp), dimension(3) :: scaling
    dx = (2._dp * real(rng(), dp) - 1._dp) * maxscaling
    dy = (2._dp * real(rng(), dp) - 1._dp) * maxscaling
    dz = (2._dp * real(rng(), dp) - 1._dp) * maxscaling
    scaling = (/(getx(simbox) + dx) / getx(simbox), &
    (gety(simbox) + dy) / gety(simbox), &
    (getz(simbox) + dz) / getz(simbox)/)    
    call setx(simbox, getx(simbox) + dx)  
    call sety(simbox, gety(simbox) + dy)  
    call setz(simbox, getz(simbox) + dz)  
  end function

  subroutine box_writetostdio(simbox)
    type(boxdat), intent(in) :: simbox
    write(*,*) simbox
  end subroutine

end module box
