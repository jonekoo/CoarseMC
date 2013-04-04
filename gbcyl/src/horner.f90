!! Calculates the value of a polynomial of order n-1 in the point x. 
!! Coefficients of the polynomial are given in table a, where 
!! a(i) is the coefficient of the x^(n-i) term
!! Calculation is done with the Horner method. 
pure function horner(a,n,x) 
  use nrtype
  implicit none
  real(dp), dimension(:),intent(in) :: a
  integer, intent(in) :: n
  real(dp), intent(in) :: x 
  real(dp) :: horner
  integer :: i
  horner=a(1)
  do i=2,n
    horner=horner*x+a(i)
  end do 
  return
end function horner
