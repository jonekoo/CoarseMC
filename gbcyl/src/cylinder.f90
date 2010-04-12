module cylinder
use nrtype
use box
implicit none
  
type cylinderdat
  private
  type(boxdat) :: b
end type cylinderdat

interface minimage
  module procedure cylinder_minimage
end interface

interface mindistance
  module procedure cylinder_mindistance
end interface

interface scale
  module procedure cyl_scale
end interface

interface writetostdio
  module procedure cylinder_writetostdio
end interface

interface setx
  module procedure cylinder_setx
end interface

interface sety
  module procedure cylinder_sety
end interface

interface setz
  module procedure cylinder_setz
end interface

interface getx
  module procedure cylinder_getx
end interface

interface gety
  module procedure cylinder_gety
end interface

interface getz
  module procedure cylinder_getz
end interface

interface isxperiodic
  module procedure cyl_isxperiodic
end interface

interface isyperiodic
  module procedure cyl_isxperiodic
end interface

interface iszperiodic
  module procedure cyl_iszperiodic
end interface

interface makebox
  module procedure cylinder_makebox, cylinder_copy
end interface

interface volume
  module procedure cylinder_volume
end interface

interface createbox
  module procedure cylinder_createbox
end interface

interface new_cylinder
  module procedure new_cylinderdh
end interface 

  contains
  
  function new_cylinderdh(diameter, height) result(nc)
    type(cylinderdat) :: nc
    real(dp), intent(in) :: diameter, height
    call makebox(nc, diameter, height)
  end function
   
  subroutine cylinder_makebox(cylbox, diameter, height) 
    type(cylinderdat), intent(out) :: cylbox
    real(dp), intent(in) :: diameter
    real(dp), intent(in) :: height
    cylbox%b = new_box(diameter)
    call setz(cylbox, height)
    call setx(cylbox, diameter)
    call setxperiodicity(cylbox%b, .false.)
    call setyperiodicity(cylbox%b, .false.)
  end subroutine

  subroutine cylinder_copy(copy, cylbox) 
    type(cylinderdat), intent(out) :: copy
    type(cylinderdat), intent(in) :: cylbox
    copy = cylbox
  end subroutine

  subroutine cylinder_createbox(cylbox, boxstring)
    type(cylinderdat), intent(out) :: cylbox
    character(len = *), intent(in) :: boxstring
    call createbox(cylbox%b, boxstring)
  end subroutine

  elemental function radius(cylbox)
    real(dp) :: radius
    type(cylinderdat), intent(in) :: cylbox
    radius = getx(cylbox) / 2._dp
  end function

  elemental function height(cylbox)
    real(dp) :: height
    type(cylinderdat), intent(in) :: cylbox
    height = getz(cylbox)
  end function

  elemental function cylinder_getx(cylbox)
    real(dp) :: cylinder_getx
    type(cylinderdat), intent(in) :: cylbox
    cylinder_getx = getx(cylbox%b)
  end function

  elemental function cylinder_gety(cylbox)
    real(dp) :: cylinder_gety
    type(cylinderdat), intent(in) :: cylbox
    cylinder_gety = gety(cylbox%b)
  end function

  elemental function cylinder_getz(cylbox)
    real(dp) :: cylinder_getz
    type(cylinderdat), intent(in) :: cylbox
    cylinder_getz = getz(cylbox%b)
  end function

  !! In a way it makes no sense to ask cylinder of its periodicity in x 
  !! or y direction. That indicates that the objects needing this 
  !! routine should be written in a way that they don't need the 
  !! information or it is given to them and they don't ask for it. 
  !! 
  elemental function cyl_isxperiodic(cylbox) result(isperiodic)
    type(cylinderdat), intent(in) :: cylbox
    logical :: isperiodic
    isperiodic = .false.
  end function

  elemental function cyl_iszperiodic(cylbox) result(isperiodic)
    type(cylinderdat), intent(in) :: cylbox
    logical :: isperiodic
    isperiodic = iszperiodic(cylbox%b)
  end function

  pure subroutine cylinder_setx(cylbox, x)
    type(cylinderdat), intent(inout) :: cylbox
    real(dp), intent(in) :: x
    call setx(cylbox%b, x)
    call sety(cylbox%b, x)
  end subroutine

  pure subroutine cylinder_sety(cylbox, y)
    type(cylinderdat), intent(inout) :: cylbox
    real(dp), intent(in) :: y
    call setx(cylbox%b, y)
    call sety(cylbox%b, y)
  end subroutine

  pure subroutine cylinder_setz(cylbox, z)
    type(cylinderdat), intent(inout) :: cylbox
    real(dp), intent(in) :: z
    call setz(cylbox%b, z)
  end subroutine

  pure function cylinder_volume(cylbox)
    intrinsic atan
    real(dp) :: cylinder_volume
    type(cylinderdat), intent(in) :: cylbox
    real(dp) :: pi
    pi = 4._dp * atan(1._dp)
    cylinder_volume = pi * radius(cylbox)**2 * height(cylbox)
  end function

  pure function cylinder_minimage(cylbox, r1, r2) result(mini)
    real(dp), dimension(3) :: mini
    type(cylinderdat), intent(in) :: cylbox
    real(dp), dimension(3), intent(in) :: r1
    real(dp), dimension(3), intent(in) :: r2
    mini = minimage(cylbox%b, r1, r2)
  end function

  pure function cylinder_mindistance(cylbox, r1, r2) result(mind)
    real(dp) :: mind
    type(cylinderdat), intent(in) :: cylbox
    real(dp), dimension(3), intent(in) :: r1
    real(dp), dimension(3), intent(in) :: r2
    mind = mindistance(cylbox%b, r1, r2)
  end function

  function cyl_scale(cylbox, maxscaling, rng) result(scaling)
    type(cylinderdat), intent(inout) :: cylbox
    real(dp), intent(in) :: maxscaling
    real(dp) :: dz
    interface 
      function rng()
        real(8) :: rng
      end function rng
    end interface
    real(dp), dimension(3) :: scaling
    dz = (2._dp * real(rng(), dp) - 1._dp) * maxscaling
    scaling = (/1._dp, 1._dp, (getz(cylbox) + dz) / getz(cylbox)/)  
    call setz(cylbox, getz(cylbox) + dz)
  end function

  subroutine cylinder_writetostdio(cylbox)
    type(cylinderdat) :: cylbox
    call writetostdio(cylbox%b)
  end subroutine

end module cylinder
