module jacobiwrapper
!! This module's purpose is to wrap the single precision jacobi routine to a
!! double precision interface. With this there need not to be explicit 
!! conversions in each module using the jacobi routine. Also the nr module and
!! the jacobi routine itself remain intact. 
!!
use num_kind
use nr, only: jacobi_sp => jacobi
implicit none

INTERFACE jacobi
  module procedure jacobi_d
END INTERFACE

contains

!! @see documentation in the book Numerical Recipes in Fortran 90. 
!! This subroutine wraps the single precision routine in the nr module.
!!
SUBROUTINE jacobi_d(a,d,v,nrot)
  INTEGER(I4B), INTENT(OUT) :: nrot
  REAL(DP), DIMENSION(:), INTENT(OUT) :: d
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: v
  REAL(SP), DIMENSION(3) :: d_sp
  REAL(SP), DIMENSION(3,3) :: a_sp
  REAL(SP), DIMENSION(3,3) :: v_sp
  a_sp = real(a, sp)
  call jacobi_sp(a_sp, d_sp, v_sp, nrot)
  a = real(a_sp, dp)
  d = real(d_sp, dp)
  v = real(v_sp, dp)
END SUBROUTINE

end module
