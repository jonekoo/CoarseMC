module orientational_ordering
  use nrtype
  use particle, only: particledat, orientation
  use jacobiwrapper, only: jacobi
  implicit none

  public :: eigens
  public :: orientation_parameter
  public :: orientation_parameter_v
  public :: orientation_tensor
  public :: cycle_largest_to_3rd

  PRIVATE

  interface eigens
    module procedure eigens_dv
  end interface eigens

  interface orientation_parameter
    module procedure orientation_parameter_s
  end interface orientation_parameter
  
contains 

!! Returns the orientation parameter of @p particles with respect to 
!! @p vector.
pure function orientation_parameter_v(particles, vector) result(res)
  type(particledat), intent(in) :: particles(:)
  real(dp), intent(in) :: vector(3)
  real(dp) :: res
  integer :: i
  res = 0._dp
  do i = 1, size(particles)
     res = res + dot_product(orientation(particles(i)), vector)**2
  end do
  res = 1.5_dp * res/size(particles) - 0.5_dp
end function orientation_parameter_v


!! Returns the largest eigenvalue of the orientational ordering tensor and
!! the eigenvector corresponding to that value. 
!!
!! pre-conditions: 
!! 1. @p particles has @p n_particles > 0 particles.
!! 
!! @p particles the particles for which the orientational ordering tensor 
!! is calculated.
!! @p n_particles the number of particles.
!! @p value to be assigned as the orientation parameter.
!! @p director the eigenvector corresponding to @p value. 
!!
subroutine orientation_parameter_s(particles, n_particles, value, director)
  implicit none
  intrinsic maxloc
  intrinsic maxval
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: n_particles
  real(dp), intent(out) :: value
  real(dp), dimension(3), intent(out) :: director
  integer, dimension(1) :: max_value_position
  integer, parameter :: n_dimensions = 3
  real(dp), dimension(3) :: values
  real(dp), dimension(3, 3) :: vectors
  call eigens_dv(particles, n_particles, values, vectors)
  value = maxval(values)
  max_value_position = maxloc(values)
  director(1:n_dimensions) = vectors(1:n_dimensions, max_value_position(1))
end subroutine orientation_parameter_s



subroutine eigens_dv(particles, n_particles, values, vectors)
  implicit none
  type(particledat), dimension(:) :: particles
  integer, intent(in) :: n_particles
  integer :: nrot
  real(dp),dimension(3,3), intent(out) :: vectors
  real(dp),dimension(3), intent(out) :: values
  real(dp),dimension(3,3) :: tensor
  call orientation_tensor(particles, n_particles, tensor)
  call jacobi(tensor, values, vectors, nrot)    
end subroutine eigens_dv


!! Makes a cyclic permutation to values and vectors so that the largest of
!! values is values(3) and the corresponding vector is vectors(:, 3)
!! 
subroutine cycle_largest_to_3rd(values, vectors)
  real(dp), intent(inout) :: values(3)
  real(dp), intent(inout) :: vectors(3, 3)
  real(dp) :: temp
  real(dp) :: temp_vector(3)
  integer :: j, i, offset, max_value_position(1)
  max_value_position = maxloc(values)
  offset = 3 - max_value_position(1)
  if (offset /= 0) then
     do i = 1, offset
        !! One cyclic permutation:
        temp = values(3)
        temp_vector = vectors(:, 3)
        do j = 2, 1, -1
           values(j + 1) = values(j)
           vectors(:, j + 1) = vectors(:, j)
        end do
        values(1) = temp
        vectors(:, 1) = temp_vector
     end do
  end if
end subroutine cycle_largest_to_3rd


!! Sorts the matrix columns by corresponding value in increasing order.
!! Implements an insertion sort as presented in introduction to algorithms by
!! T. H. Cormen et al. 
!! 
subroutine sort_by_value(values, vectors)
  real(dp), intent(inout) :: values(3)
  real(dp), intent(inout) :: vectors(3, 3)
  integer :: j, i
  real(dp) :: key
  real(dp) :: data(3)
  do j = 2, 3
     key = values(j)
     data = vectors(:, j)
     i = j - 1 
     do while(i > 0 .and. values(i) > key)
        values(i + 1) = values(i)
        vectors(:, i + 1) = vectors(:, i)
        i = i - 1
     end do
     values(i + 1) = key
     vectors(:, i + 1) = data
  end do
end subroutine sort_by_value


pure subroutine orientation_tensor(particles, n_particles, tensor)
  implicit none
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) ::  n_particles
  real(dp), dimension(3, 3), intent(out) :: tensor
  real(dp) :: Sxx, Sxy, Sxz
  real(dp) :: Syx, Syy, Syz
  real(dp) :: Szx, Szy, Szz
  integer :: i, n_rods 
  Sxx = 0._dp; Sxy = 0._dp; Sxz = 0._dp;
  Syx = 0._dp; Syy = 0._dp; Syz = 0._dp;
  Szx = 0._dp; Szy = 0._dp; Szz = 0._dp;
  n_rods = 0;
  do i = 1, n_particles
     if (particles(i)%rod) then
        Sxx = Sxx + 3._dp*particles(i)%ux*particles(i)%ux - 1._dp;
        Sxy = Sxy + 3._dp*particles(i)%ux*particles(i)%uy;
        Sxz = Sxz + 3._dp*particles(i)%ux*particles(i)%uz;
        
        Syx = Syx + 3._dp*particles(i)%uy*particles(i)%ux;
        Syy = Syy + 3._dp*particles(i)%uy*particles(i)%uy - 1._dp;
        Syz = Syz + 3._dp*particles(i)%uy*particles(i)%uz;
        
        Szx = Szx + 3._dp*particles(i)%uz*particles(i)%ux;
        Szy = Szy + 3._dp*particles(i)%uz*particles(i)%uy;
        Szz = Szz + 3._dp*particles(i)%uz*particles(i)%uz - 1._dp;
        
        n_rods = n_rods+1;
     end if
  end do
    
  Sxx=Sxx/(2._dp*real(n_rods, dp)); Sxy=Sxy/(2._dp*real(n_rods, dp)); Sxz=Sxz/(2._dp*real(n_rods, dp));
  Syx=Syx/(2._dp*real(n_rods, dp)); Syy=Syy/(2._dp*real(n_rods, dp)); Syz=Syz/(2._dp*real(n_rods, dp));
  Szx=Szx/(2._dp*real(n_rods, dp)); Szy=Szy/(2._dp*real(n_rods, dp)); Szz=Szz/(2._dp*real(n_rods, dp));
  
  tensor(1,1)=Sxx; tensor(1,2)=Sxy; tensor(1,3)=Sxz;
  tensor(2,1)=Syx; tensor(2,2)=Syy; tensor(2,3)=Syz;
  tensor(3,1)=Szx; tensor(3,2)=Szy; tensor(3,3)=Szz;
end subroutine orientation_tensor

end module
