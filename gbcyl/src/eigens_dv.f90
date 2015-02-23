subroutine eigensystem(vectors, values)
  use nrtype, only: dp
  use jacobiwrapper, only: jacobi
  implicit none
  real(dp),dimension(3, 3), intent(inout) :: vectors
  real(dp),dimension(3), intent(out) :: values
  real(dp),dimension(3, 3) :: tensor
  integer :: nrot
  tensor = vectors
  call jacobi(tensor, values, vectors, nrot)    
end subroutine


