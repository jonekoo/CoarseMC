FUNCTION locate(xx,x)
  USE nrtype
  IMPLICIT NONE
  real(sp), dimension(:), intent(in) :: xx
  real(sp), intent(in) :: x
  integer(I4B) :: locate
  !! Given an array xx(1:N), and given a value x, returns a value j such that
  !! x is between xx(j) and xx(j+1). xx must be monotonic, either increasing
  !! or decreasing. j=0 or j=N is returned to indicate that x is out of range.
  integer(I4B) :: n,jl,jm,ju
  logical :: ascnd
  n=size(xx)
  ascnd=(xx(n) >= xx(1))   !! True if ascending order of table, false otherwise
  jl=0                     !! Initialize lower
  ju=n+1                   !! and upper limits.
  do 
    if (ju-jl <= 1) exit   !! Repeat until this condition is satisfied.
    jm=(ju+jl)/2           !! Compute a midpoint
    if (ascnd .eqv. (x >=xx(jm))) then
      jl=jm                !! and replace either the lower limit
    else 
      ju=jm                !! or the upper limit, as appropriate.
    end if
  end do
  if (x==xx(1)) then       !! Then set the output, being careful with the 
    locate=1               !! endpoints.
  else if (x==xx(n)) then
    locate=n-1
  else 
    locate=jl
  end if
END FUNCTION locate
  
   
