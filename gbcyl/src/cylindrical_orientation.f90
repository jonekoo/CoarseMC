module cylindrical_orientation
  use num_kind, only: dp
  use m_particledat
  use utils, only: unitvec
  use orientational_ordering, only: orientational_ordering_eigens => eigens
  implicit none  

  contains
  subroutine eigens(particles, n_particles, values, vectors)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp), dimension(3), intent(out) :: values
    real(dp), dimension(3, 3), intent(out) :: vectors
    type(particledat), dimension(:), allocatable :: temp_particles
    integer :: i
    allocate(temp_particles(n_particles))
    temp_particles(1:n_particles) = particles(1:n_particles)
    do i = 1, n_particles
      if(temp_particles(i)%rod) then 
         call setorientation(temp_particles(i), &
              unitvec(temp_particles(i)%orientation(), &
             temp_particles(i)%position()))
      end if
    end do
    call orientational_ordering_eigens(temp_particles, n_particles, &
      & values, vectors)
    if (allocated(temp_particles)) deallocate(temp_particles)
  end subroutine
end module 
