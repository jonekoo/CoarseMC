module tau1_negative
use nrtype
implicit none

public :: tau1_n
public :: tau1
public :: init
private


interface tau1
  module procedure tau1_r
end interface

!!interface 
!!  function tau1(rs, n, d)
!!    use nrtype, only: sp
!!    real(sp), dimension(:), intent(in) :: rs
!!    integer, intent(in) :: n
!!    real(sp), intent(in) :: d
!!  end function
!!end interface

real(sp), dimension(:), allocatable, save :: rs_
integer, save :: n_particles_ = 0

contains

subroutine init(rs, n_particles)
  real(sp), dimension(:), intent(in) :: rs
  integer, intent(in) :: n_particles
  integer :: allocstat
  allocstat = 0
  if(n_particles /= n_particles_) then
    n_particles_ = n_particles
    if(allocated(rs_)) deallocate(rs_)
  end if
  if(.not. allocated(rs_)) then
    allocate(rs_(n_particles_), stat = allocstat)
  end if
  if(allocstat /= 0) stop 'Allocation failed in module tau1_negative'
  rs_(1:n_particles_) = rs(1:n_particles_)
end subroutine

function tau1_n(layer_spacing) 
  implicit none
  real(sp), intent(in) :: layer_spacing
  real(sp) :: tau1_n
  tau1_n = -tau1(rs_, n_particles_, layer_spacing)   
end function

function tau1_r(rs, n, d) result(t1)
  real(sp), dimension(:), intent(in) :: rs
  integer, intent(in) :: n
  real(sp), intent(in) :: d
  integer :: k 
  real(sp) :: pi
  real(sp) :: t1
  complex(spc) :: t1c
  complex(spc), parameter :: i=(0._sp, 1._sp)
  pi = 4._sp * atan(1._sp)
  t1c = (0._sp, 0._sp)
  do k = 1, n
    t1c = t1c + exp(i * cmplx(2._sp * pi * rs(k) / d, 0._sp, spc))
  end do
  t1c = t1c / cmplx(n, 0._sp, spc)
  t1 = abs(t1c)
end function

end module
