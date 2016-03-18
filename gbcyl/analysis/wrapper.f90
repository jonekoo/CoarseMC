subroutine wrap_particledat(arr, n, p2, director)
  use m_particledat, only: particledat
  use orientational_ordering, only: orientation_parameter
  use num_kind, only: dp
  use utils, only: fmt_char_dp
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: arr(6, n)
  !f2py real(8) intent(in), depend(n) :: arr
  type(particledat), allocatable :: particles(:)
  real(8), intent(out) :: p2
  real(8), intent(out) :: director(3)
  integer :: i
  allocate(particles(size(arr, 2)))
  do i = 1, size(arr, 2)
     particles(i)%x=arr(1, i)
     particles(i)%y=arr(2, i)
     particles(i)%z=arr(3, i)
     particles(i)%ux=arr(4, i)
     particles(i)%uy=arr(5, i)
     particles(i)%uz=arr(6, i)
  end do
  particles%rod = .true.
  call orientation_parameter(particles, size(particles), p2, director)
end subroutine wrap_particledat


