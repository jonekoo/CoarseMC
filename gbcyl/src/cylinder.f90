module cylinder
use nrtype
use box
implicit none
  
type cylinderdat
  private
  type(boxdat) :: b
end type cylinderdat

interface min_image
  module procedure cylinder_min_image
end interface

interface min_distance
  module procedure cylinder_min_distance
end interface

interface scale
  module procedure cyl_scale
end interface

interface write_to_stdio
  module procedure cylinder_write_to_stdio
end interface

interface set_x
  module procedure cylinder_set_x
end interface

interface set_y
  module procedure cylinder_set_y
end interface

interface set_z
  module procedure cylinder_set_z
end interface

interface get_x
  module procedure cylinder_get_x
end interface

interface get_y
  module procedure cylinder_get_y
end interface

interface get_z
  module procedure cylinder_get_z
end interface

interface is_xperiodic
  module procedure cyl_is_xperiodic
end interface

interface is_yperiodic
  module procedure cyl_is_xperiodic
end interface

interface is_zperiodic
  module procedure cyl_is_zperiodic
end interface

interface make_box
  module procedure cylinder_make_box, cylinder_copy
end interface

interface volume
  module procedure cylinder_volume
end interface

interface create_box
  module procedure cylinder_create_box
end interface

interface new_cylinder
  module procedure new_cylinder_dh
end interface 

  contains
  
  function new_cylinder_dh(diameter, height) result(nc)
    type(cylinderdat) :: nc
    real(dp), intent(in) :: diameter, height
    call make_box(nc, diameter, height)
  end function
   
  subroutine cylinder_make_box(cylbox, diameter, height) 
    type(cylinderdat), intent(out) :: cylbox
    real(dp), intent(in) :: diameter
    real(dp), intent(in) :: height
    cylbox%b = new_box(diameter)
    call set_z(cylbox, height)
    call set_x(cylbox, diameter)
    call set_xperiodicity(cylbox%b, .false.)
    call set_yperiodicity(cylbox%b, .false.)
  end subroutine

  subroutine cylinder_copy(copy, cylbox) 
    type(cylinderdat), intent(out) :: copy
    type(cylinderdat), intent(in) :: cylbox
    copy = cylbox
  end subroutine

  subroutine cylinder_create_box(cylbox, box_string)
    type(cylinderdat), intent(out) :: cylbox
    character(len = *), intent(in) :: box_string
    call create_box(cylbox%b, box_string)
  end subroutine

  elemental function radius(cylbox)
    real(dp) :: radius
    type(cylinderdat), intent(in) :: cylbox
    radius = get_x(cylbox) / 2._dp
  end function

  elemental function height(cylbox)
    real(dp) :: height
    type(cylinderdat), intent(in) :: cylbox
    height = get_z(cylbox)
  end function

  elemental function cylinder_get_x(cylbox)
    real(dp) :: cylinder_get_x
    type(cylinderdat), intent(in) :: cylbox
    cylinder_get_x = get_x(cylbox%b)
  end function

  elemental function cylinder_get_y(cylbox)
    real(dp) :: cylinder_get_y
    type(cylinderdat), intent(in) :: cylbox
    cylinder_get_y = get_y(cylbox%b)
  end function

  elemental function cylinder_get_z(cylbox)
    real(dp) :: cylinder_get_z
    type(cylinderdat), intent(in) :: cylbox
    cylinder_get_z = get_z(cylbox%b)
  end function

  !! In a way it makes no sense to ask cylinder of its periodicity in x 
  !! or y direction. That indicates that the objects needing this 
  !! routine should be written in a way that they don't need the 
  !! information or it is given to them and they don't ask for it. 
  !! 
  elemental function cyl_is_xperiodic(cylbox) result(is_periodic)
    type(cylinderdat), intent(in) :: cylbox
    logical :: is_periodic
    is_periodic = .false.
  end function

  elemental function cyl_is_zperiodic(cylbox) result(is_periodic)
    type(cylinderdat), intent(in) :: cylbox
    logical :: is_periodic
    is_periodic = is_zperiodic(cylbox%b)
  end function

  pure subroutine cylinder_set_x(cylbox, x)
    type(cylinderdat), intent(inout) :: cylbox
    real(dp), intent(in) :: x
    call set_x(cylbox%b, x)
    call set_y(cylbox%b, x)
  end subroutine

  pure subroutine cylinder_set_y(cylbox, y)
    type(cylinderdat), intent(inout) :: cylbox
    real(dp), intent(in) :: y
    call set_x(cylbox%b, y)
    call set_y(cylbox%b, y)
  end subroutine

  pure subroutine cylinder_set_z(cylbox, z)
    type(cylinderdat), intent(inout) :: cylbox
    real(dp), intent(in) :: z
    call set_z(cylbox%b, z)
  end subroutine

  pure function cylinder_volume(cylbox)
    intrinsic atan
    real(dp) :: cylinder_volume
    type(cylinderdat), intent(in) :: cylbox
    real(dp) :: pi
    pi = 4._dp * atan(1._dp)
    cylinder_volume = pi * radius(cylbox)**2 * height(cylbox)
  end function

  pure function cylinder_min_image(cylbox, r1, r2) result(min_i)
    real(dp), dimension(3) :: min_i
    type(cylinderdat), intent(in) :: cylbox
    real(dp), dimension(3), intent(in) :: r1
    real(dp), dimension(3), intent(in) :: r2
    min_i = min_image(cylbox%b, r1, r2)
  end function

  pure function cylinder_min_distance(cylbox, r1, r2) result(min_d)
    real(dp) :: min_d
    type(cylinderdat), intent(in) :: cylbox
    real(dp), dimension(3), intent(in) :: r1
    real(dp), dimension(3), intent(in) :: r2
    min_d = min_distance(cylbox%b, r1, r2)
  end function

  function cyl_scale(cylbox, max_scaling, rng) result(scaling)
    type(cylinderdat), intent(inout) :: cylbox
    real(dp), intent(in) :: max_scaling
    real(dp) :: dz
    interface 
      function rng()
        real(8) :: rng
      end function rng
    end interface
    real(dp), dimension(3) :: scaling
    dz = (2._dp * real(rng(), dp) - 1._dp) * max_scaling
    scaling = (/1._dp, 1._dp, (get_z(cylbox) + dz) / get_z(cylbox)/)  
    call set_z(cylbox, get_z(cylbox) + dz)
  end function

  subroutine cylinder_write_to_stdio(cylbox)
    type(cylinderdat) :: cylbox
    call write_to_stdio(cylbox%b)
  end subroutine

end module cylinder
