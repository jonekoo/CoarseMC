module tau1_negative
use num_kind
implicit none

public :: tau1_n
public :: tau1
public :: tau1_c
public :: init
private


interface tau1
  module procedure tau1_r
end interface

real(sp), dimension(:), allocatable, save :: rs_


contains


!> Initializes the module with @p rs as the particle positions 
!! projected to the reference direction (usually the liquid crystal
!! director).
!!
subroutine init(rs)
  real(sp), dimension(:), intent(in) :: rs
  integer :: allocstat
  allocstat = 0
  if (allocated(rs_)) then
     if(size(rs) /= size(rs_)) deallocate(rs_)
  end if
  if(.not. allocated(rs_)) then
    allocate(rs_(size(rs)), stat = allocstat)
  end if
  if(allocstat /= 0) stop 'Allocation failed in module tau1_negative'
  rs_(:) = rs(:)
end subroutine

!> Returns -tau_1 with a given @p layer_spacing. Meant to be used with a
!! minimization routine to find the value of the order parameter which
!! is the maximum value of this.
!!
!! @param layer_spacing the distance between layers.
!!
!! @return the value of -tau_1.
!!
function tau1_n(layer_spacing) 
  implicit none
  real(sp), intent(in) :: layer_spacing
  real(sp) :: tau1_n
  tau1_n = -tau1(rs_, layer_spacing)   
end function

!> Returns the absolute value of tau_1 with projected particle positions
!! @p rs and layer distance @p d.
!!
function tau1_r(rs, d) result(t1)
  real(sp), dimension(:), intent(in) :: rs
  real(sp), intent(in) :: d
  real(sp) :: t1
  t1 = abs(tau1_c(rs, d))
end function


!> Returns the complex value of tau_1 with given projected particle
!! positions @p rs and layer distance @p d.
!!
function tau1_c(rs, d) result(t1c)
  real(sp), dimension(:), intent(in) :: rs
  real(sp), intent(in) :: d
  integer :: k 
  real(sp) :: pi
  complex(spc) :: t1c
  complex(spc), parameter :: i = (0, 1)
  pi = 4._sp * atan(1._sp)
  t1c = (0._sp, 0._sp)
  do k = 1, size(rs)
    t1c = t1c + exp(i * 2._sp * pi * rs(k) / d)
  end do
  t1c = t1c / size(rs)
end function

end module
