!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Some prototyping below here.                                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Maps @p func to all @p particles. The cell list @p sl is used to do
!! this in parallel in a checkerboard order. 
!! 
!! @param particles is the array of particles to move.
!! @param sl the cell list presenting the decomposition.
!! 
subroutine map(sl, particles, routine)
  implicit none
  type(particledat), intent(inout) :: particles(:)
  type(simplelist), intent(in) :: sl
  interface
     !! Voisiko olla elemental?
     subroutine routine(p)
       use particle
       type(particledat), intent(inout) :: p
     end subroutine routine
  end interface
       
  !$ integer :: n_threads = 1
  integer :: thread_id = 0
  !! You may be tempted to make the allocatable arrays automatic, but there's
  !! no performance gained and depending on the system large automatic arrays
  !! may cause a stack overflow.
  type(particledat), allocatable :: temp_particles(:)
  integer :: n_cell ! particles in the cell where particles are modified.
  integer :: j, ix, iy, iz, jx, jy, jz
  !! Loop over cells. This can be thought of as looping through a 
  !! 2 x 2 x 2 cube of cells.
  !$OMP PARALLEL shared(particles, sl)& 
  !$OMP& private(thread_id, n_threads, temp_particles)
  !$ thread_id = omp_get_thread_num()
  allocate(temp_particles(maxval(sl%counts)))
  do iz=0, min(1, sl%nz-1)
     do iy=0, min(1, sl%ny-1)
        do ix=0, min(1, sl%nx-1)
           !$OMP DO collapse(3) private(j, n_cell)&
           !$OMP& schedule(dynamic)
           do jz = iz, sl%nz - 1, 2
              do jy = iy, sl%ny - 1, 2
                 do jx = ix, sl%nx - 1, 2
                    n_cell = sl%counts(jx, jy, jz)
                    temp_particles(1:n_cell) = &
                         particles(sl%indices(1:n_cell, jx, jy, jz))
                    do j = 1, n_cell
                       call routine(temp_particles(j))
                    end do
                    particles(sl%indices(1:n_cell, jx, jy, jz)) = &
                         temp_particles(1:n_cell)
                 end do
              end do
           end do
           !$OMP END DO 
           !! The end of parallelized loop forces an implicit barrier.
           !! Memory view of the threads is also synchronized here.
        end do
        !$OMP BARRIER
     end do
     !$OMP BARRIER
  end do
  if (allocated(temp_particles)) deallocate(temp_particles)
  !$OMP END PARALLEL
end subroutine map




