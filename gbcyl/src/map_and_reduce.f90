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




!> Computes the total @p energy of the system. Interactions are computed
!! cell-by-cell for the @p particles in the cell list @p sl. @p simbox
!! is the simulation cell. If any two particles are too close to each
!! other, @p overlap is true. 
subroutine pair_reduce(sl, particles, accumulate, energy, overlap)
  type(simplelist), intent(in) :: sl
  !type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: energy
  interface
     subroutine accumulate(p, others, energy, overlap)
       use particle
       implicit none
       type(particledat), intent(in) :: p, others(:)
       real(dp), intent(out) :: energy
       logical, intent(out) :: overlap
     end subroutine accumulate
  end interface
       
  logical, intent(out) :: overlap 
  integer :: i, j, ix, iy, iz
  logical :: mask(size(particles))
  real(dp) :: energy_j
  logical :: overlap_j
  integer :: helper(size(particles))
  integer, allocatable :: temp_helper(:)
  integer :: n_mask
  integer :: temp_j
  type(particledat), allocatable :: temp_particles(:)

  helper = (/(i, i=1, size(particles))/) !! ifort vectorizes
  energy = 0._dp
  overlap = .false.
  
  !$OMP PARALLEL default(shared) reduction(+:energy) reduction(.or.:overlap)& 
  !$OMP& private(energy_j, overlap_j, i, j, mask, temp_particles, temp_j, temp_helper, n_mask)
  allocate(temp_particles(size(particles)), temp_helper(size(particles))) 
  !$OMP DO collapse(3) schedule(dynamic)
  do ix=0, sl%nx-1 
  do iy=0, sl%ny-1
  do iz=0, sl%nz-1
    call simplelist_nbrmask(sl, simbox, ix, iy, iz, mask)
    n_mask = count(mask) 
    temp_particles(1:n_mask) = pack(particles, mask)
    temp_helper(1:n_mask) = pack(helper, mask)
    do i = 1, sl%counts(ix, iy, iz)
      j = sl%indices(i, ix, iy, iz)
      ! Find position of particles(j) in temp_particles:
      do temp_j = 1, n_mask 
         if(temp_helper(temp_j) == j) exit
      end do
      call accumulate(temp_particles(temp_j), &
           temp_particles(temp_j + 1:n_mask), energy_j, overlap_j)
      overlap = overlap .or. overlap_j 
      energy = energy + energy_j
    end do 
  end do
  end do
  end do
  !$OMP END DO 
  if (allocated(temp_particles)) deallocate(temp_particles)
  if (allocated(temp_helper)) deallocate(temp_helper)
  !$OMP END PARALLEL  
end subroutine pair_reduce


