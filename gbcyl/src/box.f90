module box
use nrtype
implicit none

type boxdat
  real(dp) :: lx
  real(dp) :: ly
  real(dp) :: lz
  logical :: xperiodic
  logical :: yperiodic
  logical :: zperiodic    
end type boxdat

contains 

  !! :TODO: Assignment operator, if needed. 


  !! Returns a "well-defined" periodic cubic box with sidelengths @p side.
  !! This routine is meant to be a initializing routine after which 
  !! customizations can be done.
  !!
  !! @p side the length of one side.
  !!
  !! @return a cubic box.
  !!  
  type(boxdat) function make_box(side)
  implicit none
  real(dp), intent(in) :: side
    make_box%lx = side
    make_box%ly = side
    make_box%lz = side
    make_box%xperiodic = .true.
    make_box%yperiodic = .true.
    make_box%zperiodic = .true.
  end function



  real(dp) function box_x(a_box)
    type(boxdat), intent(in) :: a_box
    box_x = a_box%lx
  end function



  real(dp) function box_y(a_box)
    type(boxdat), intent(in) :: a_box
    box_y = a_box%ly
  end function



  real(dp) function box_z(a_box)
    type(boxdat), intent(in) :: a_box
    box_z = a_box%lz
  end function



  !!subroutine box_x(a_box, lx)
  !!  type(boxdat, intent(in) :: a_box
  !!  a_box%lx = x
  !!end subroutine box_x

end module box
