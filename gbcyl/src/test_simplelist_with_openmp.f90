program simplelist_with_openmp
use class_simplelist
!$ use omp_lib
implicit none
type(simplelist) :: sl
integer, parameter :: nx = 4, ny = 4, nz = 4
integer, parameter :: n = nx * ny * nz * 2
integer :: i
integer :: thread_id = 0
integer :: indices(n, max(nx/2,1) * max(ny/2, 1) * max(nz/2, 1))
integer :: ix, iy, iz
sl%nx = nx
sl%ny = ny
sl%nz = nz
call simplelist_allocate(sl, n)
sl%indices = 0
sl%indices(1:2, :, :, :) = reshape((/(i, i=1, n)/), (/2, nx, ny, nz/)) 
!$OMP PARALLEL
!$ thread_id = omp_get_thread_num()
!$ if (thread_id == 0) write(*, *) "Number of threads: ", omp_get_num_threads()
!$OMP DO
do i = 0, sl%nz - 1, 2
  write(*, *) "thread:", thread_id, "indices:", sl%indices(:, 0, 0, i)
end do
!$OMP END DO

!! Loop over cells. This can be thought of as looping through a 2 x 2 x 2 cube
!! of cells.
do ix=0, min(1, sl%nx-1)
do iy=0, min(1, sl%ny-1)
do iz=0, min(1, sl%nz-1)
indices = reshape(sl%indices(:, ix:sl%nx-1:2, iy:sl%ny-1:2, iz:sl%nz-1:2), (/n, size(indices(1,:))/))
!$ thread_id = omp_get_thread_num()
!$OMP DO
do i = 1, size(indices(1,:))
  write(*, *) "thread:", thread_id, ix, iy, iz, pack(indices(:,i), indices(:,i) > 0)
end do
!$OMP END DO
end do 
end do
end do
!$OMP END PARALLEL

end program
