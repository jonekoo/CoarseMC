program test_dsyev
implicit none
integer, parameter :: n = 3
integer, parameter :: lwork = max(1, 3 * n - 1)
real(8) :: a(n, n), w(n), work(lwork)
integer :: info
interface 
   subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
     character, intent(in) :: jobz
     !! = 'N':  Compute eigenvalues only;
     !! = 'V':  Compute eigenvalues and eigenvectors.
     character, intent(in) :: uplo
     !! = 'U':  Upper triangle of A is stored;
     !! = 'L':  Lower triangle of A is stored.
     integer, intent(in) :: n
     !! The order of the matrix A.  N >= 0.
     integer, intent(in) :: lda
     !! The leading dimension of the array A.  LDA >= max(1,N).
     double precision, dimension(lda, *), intent(inout) :: a
     !! A is DOUBLE PRECISION array, dimension (LDA, N)
     !! On entry, the symmetric matrix A.  If UPLO = 'U', the
     !! leading N-by-N upper triangular part of A contains the
     !! upper triangular part of the matrix A.  If UPLO = 'L',
     !! the leading N-by-N lower triangular part of A contains
     !! the lower triangular part of the matrix A.
     !! On exit, if JOBZ = 'V', then if INFO = 0, A contains the
     !! orthonormal eigenvectors of the matrix A.
     !! If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
     !! or the upper triangle (if UPLO='U') of A, including the
     !! diagonal, is destroyed.
     integer, intent(in) :: lwork
     !! LWORK is INTEGER
     !! The length of the array WORK.  LWORK >= max(1,3*N-1).
     !! For optimal efficiency, LWORK >= (NB+2)*N,
     !! where NB is the blocksize for DSYTRD returned by ILAENV.
     !!
     !! If LWORK = -1, then a workspace query is assumed; the routine
     !! only calculates the optimal size of the WORK array, returns
     !! this value as the first entry of the WORK array, and no error
     !! message related to LWORK is issued by XERBLA.
     double precision, dimension(*), intent(out) :: w
     !! If INFO = 0, the eigenvalues in ascending order.
     double precision, dimension(*), intent(out) :: work
     !! WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
     !! On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
     integer, intent(out) :: info
     !! = 0:  successful exit
     !! < 0:  if INFO = -i, the i-th argument had an illegal value
     !! > 0:  if INFO = i, the algorithm failed to converge; i
     !! off-diagonal elements of an intermediate tridiagonal
     !! form did not converge to zero.
   end subroutine dsyev
end interface

call init_matrix(a)
a(1, 1) = -1.
call print_matrix(a)
call dsyev('V', 'U', n, a, n, w, work, lwork, info) 
write(*, *) 'info = ', info
call print_matrix(a)


contains 

subroutine init_matrix(a)
  real(8), intent(out) :: a(:, :)
  integer :: n, i, j 
  n = size(a, 1) 
  do i = 1, n
     do j = 1, n
        a(i, j) = i * j
     end do
  end do
end subroutine


subroutine print_matrix(a)
  real(8), intent(in) :: a(:, :)
  integer :: n, i, j
  n = size(a, 1) 
  write(*, *) 'a = '
  do i = 1, n
     do j = 1, n
        write(*, fmt='((G12.6),1X)', advance='no') a(i, j)
     end do
     write(*, '(/)')
  end do
end subroutine
end program
