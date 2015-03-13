module cylindrical_layerorientation
  use nrtype, only: dp
  use particle
  use utils
  use orientational_ordering, only: orientational_ordering_eigens => eigens
  use layernormal
  use class_poly_box
  implicit none  

  real(dp), parameter :: sigma0 = 1._dp
  real(dp), save :: cutoff = 2._dp * sigma0
  
  contains

  subroutine eigens(simbox, particles, n_particles, values, vectors)
    type(poly_box), intent(in) :: simbox
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
        call setorientation(temp_particles(i), unitvec(localnormal(simbox, particles, i, cutoff), position(temp_particles(i))))
      end if
    end do
    call orientational_ordering_eigens(temp_particles, n_particles, &
      & values, vectors)
    if (allocated(temp_particles)) deallocate(temp_particles)
  end subroutine



end module 
