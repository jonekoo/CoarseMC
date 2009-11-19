module box
use nrtype
use particle
implicit none

public :: boxdat
public :: new_box
public :: make_box
public :: create_box
public :: min_image
public :: min_distance
public :: get_x, get_y, get_z
public :: set_x
public :: set_y
public :: set_z
public :: volume
public :: scale
public :: set_xperiodicity
public :: set_yperiodicity
public :: set_zperiodicity
public :: write_to_stdio

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

interface min_image
  module procedure min_image_1, min_image_2
end interface

interface min_distance
  module procedure box_min_distance
end interface

interface scale
  module procedure box_scale
end interface

interface set_x
  module procedure box_set_x
end interface

interface set_y
  module procedure box_set_y
end interface

interface set_z
  module procedure box_set_z
end interface

interface get_x
  module procedure box_get_x
end interface

interface get_y
  module procedure box_get_y
end interface

interface get_z
  module procedure box_get_z
end interface

interface set_xperiodicity
  module procedure box_set_xperiodicity
end interface

interface set_yperiodicity
  module procedure box_set_yperiodicity
end interface

interface set_zperiodicity
  module procedure box_set_zperiodicity
end interface

interface write_to_stdio
  module procedure box_write_to_stdio
end interface

interface volume
  module procedure box_volume
end interface

interface make_box
  module procedure box_make_box, box_copy
end interface

interface create_box
  module procedure box_create_box
end interface

contains 

  !! :TODO: Assignment operator, if needed. 

  function new_box(side) result(b)
    type(boxdat) :: b
    real(dp), intent(in) :: side
    call box_make_box(b, side)
  end function

  !! Returns a "well-defined" periodic cubic box with sidelengths @p side.
  !! This routine is meant to be a initializing routine after which 
  !! customizations can be done.
  !!
  !! @p side the length of one side.
  !!
  !! @return a cubic box.
  !!  
  subroutine box_make_box(simbox, side) 
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

  subroutine box_create_box(bp, box_string)
    type(boxdat), intent(inout) :: bp
    character(len = *), intent(in) :: box_string
    character(len = 5) :: type_id
    integer :: ios
    read(box_string, *, iostat = ios) type_id, bp%lx, bp%ly, bp%lz, bp%xperiodic, bp%yperiodic, bp%zperiodic
    if (ios /= 0) then
      write(6, *) 'box_create_box: failed to read box with iostat = ', ios, '. Stopping.'
      stop 
    else if('box' /= type_id) then
      write(6, *) 'box_create_box: Warning converting from ', trim(adjustl(type_id)), 'to box!'
    end if
  end subroutine

  real(dp) pure function box_get_x(a_box)
    type(boxdat), intent(in) :: a_box
    box_get_x = a_box%lx
  end function

  real(dp) pure function box_get_y(a_box)
    type(boxdat), intent(in) :: a_box
    box_get_y = a_box%ly
  end function

  real(dp) pure function box_get_z(a_box)
    type(boxdat), intent(in) :: a_box
    box_get_z = a_box%lz
  end function

  subroutine box_set_x(a_box, x)
    implicit none
    type(boxdat), intent(inout) :: a_box
    real(dp), intent(in) :: x
    a_box%lx = x
  end subroutine 
  
  subroutine box_set_y(a_box, y)
    implicit none
    type(boxdat), intent(inout) :: a_box
    real(dp), intent(in) :: y
    a_box%ly = y
  end subroutine 
  
  subroutine box_set_z(a_box, z)
    implicit none
    type(boxdat), intent(inout) :: a_box
    real(dp), intent(in) :: z
    a_box%lz = z
  end subroutine 

  subroutine box_set_xperiodicity(a_box, is_periodic)
    type(boxdat), intent(inout) :: a_box
    logical, intent(in) :: is_periodic
    a_box%xperiodic = is_periodic
  end subroutine

  subroutine box_set_yperiodicity(a_box, is_periodic)
    type(boxdat), intent(inout) :: a_box
    logical, intent(in) :: is_periodic
    a_box%yperiodic = is_periodic
  end subroutine

  subroutine box_set_zperiodicity(a_box, is_periodic)
    type(boxdat), intent(inout) :: a_box
    logical, intent(in) :: is_periodic
    a_box%zperiodic = is_periodic
  end subroutine

  pure function box_volume(simbox)
    real(dp) :: box_volume
    type(boxdat), intent(in) :: simbox
    box_volume = simbox%lx * simbox%ly * simbox%lz
  end function

  function box_min_distance(simbox, r_i, r_j)
    real(dp) :: box_min_distance
    type(boxdat), intent(in) :: simbox
    real(dp), dimension(3), intent(in) :: r_i
    real(dp), dimension(3), intent(in) :: r_j
    real(dp), dimension(3) :: r_ij
    r_ij = min_image(simbox, r_i, r_j)
    box_min_distance = sqrt(dot_product(r_ij, r_ij))
  end function 

  function min_image_2(simbox, r_i, r_j) result(r_ij)
    real(dp), dimension(3) :: r_ij
    type(boxdat), intent(in) :: simbox
    real(dp), dimension(3), intent(in) :: r_i
    real(dp), dimension(3), intent(in) :: r_j
    r_ij = r_j - r_i
    r_ij = min_image_1(simbox, r_ij)
  end function

  recursive function min_image_1(simbox, r)
    implicit none
    real(dp), dimension(3) :: min_image_1
    type(boxdat), intent(in) :: simbox
    real(dp), dimension(3), intent(in) :: r
    !! Make periodic transformations
    min_image_1 = r
    if (simbox%xperiodic) then
      min_image_1(1) = r(1) - simbox%lx * anint(r(1)/simbox%lx)
    end if
    if (simbox%yperiodic) then
      min_image_1(2) = r(2) - simbox%ly * anint(r(2)/simbox%ly)
    end if
    if (simbox%zperiodic) then
      min_image_1(3) = r(3) - simbox%lz * anint(r(3)/simbox%lz)
    end if
  end function

  subroutine box_scale(simbox, max_scaling, rng)
    type(boxdat), intent(inout) :: simbox
    real(dp), intent(in) :: max_scaling
    real(dp) :: dx, dy, dz
    interface
      function rng()
        real(8) :: rng
      end function
    end interface
    !! Make isotropic scaling. How to implement other kinds?
    dx = (2._dp * real(rng(), dp) - 1._dp) * max_scaling
    call set_x(simbox, get_x(simbox) + dx)  
    dy = (2._dp * real(rng(), dp) - 1._dp) * max_scaling
    call set_y(simbox, get_y(simbox) + dy)  
    dz = (2._dp * real(rng(), dp) - 1._dp) * max_scaling
    call set_z(simbox, get_z(simbox) + dz)  
  end subroutine

  subroutine box_write_to_stdio(simbox)
    type(boxdat), intent(in) :: simbox
    write(*,*) simbox
  end subroutine

end module box
