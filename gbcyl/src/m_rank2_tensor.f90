module m_rank2_tensor
use num_kind
implicit none

type rank2_tensor
  real(dp) :: components(3, 3)
  real(dp) :: axes(3, 3)   !! axes(1:3, 1) == local x-axis etc. 
end type

interface operator(+)
  module procedure add_components
end interface

interface operator(/)
  module procedure divide_by_int
  module procedure divide_by_real
end interface

contains


elemental function divide_by_int(tensor, int) result(t)
  type(rank2_tensor), intent(in) :: tensor
  integer, intent(in) :: int
  type(rank2_tensor) :: t
  t = tensor / real(int, dp)
end function


elemental function divide_by_real(tensor, r) result(t)
  type(rank2_tensor), intent(in) :: tensor
  real(dp), intent(in) :: r
  type(rank2_tensor) :: t
  t = tensor
  t%components = t%components / r
end function


pure function new_tensor(axes) result(t)
  real(dp), intent(in) :: axes(3, 3)
  type(rank2_tensor) :: t
  t%components=0._dp
  t%axes = axes
end function


elemental function rotate(tensor, axes) result(rotated)
  type(rank2_tensor), intent(in) :: tensor
  type(rank2_tensor), intent(in) :: axes
  type(rank2_tensor) :: rotated
  real(dp) :: cosines(3, 3)
  integer :: a, b, i, j

  do a = 1, 3
    do i = 1, 3
      cosines(a, i) = dot_product(axes%axes(1:3, a), tensor%axes(1:3, i)) 
    end do
  end do

  rotated%components = 0._dp
  do a = 1, 3
  do b = 1, 3
    do i = 1, 3
      do j = 1, 3
        rotated%components(a, b) = rotated%components(a, b) + cosines(a, i) * cosines(b, j) * tensor%components(i, j)
      end do
    end do
  end do
  end do

  rotated%axes = axes%axes
end function


function average(tensors, axes) result(res)
  type(rank2_tensor), intent(in) :: tensors(:)
  real(dp), intent(in) :: axes(3, 3)
  type(rank2_tensor) :: res
  type(rank2_tensor) :: reference
  reference%axes = axes
  res = sum(rotate(tensors, reference))
  res%components = res%components / size(tensors)
end function


function accumulate(tensors, axes) result(res)
  type(rank2_tensor), intent(in) :: tensors(:)
  real(dp), intent(in) :: axes(3, 3)
  type(rank2_tensor) :: res
  type(rank2_tensor) :: reference
  reference%axes = axes
  res = sum(rotate(tensors, reference))
end function


elemental function add_components(tensor, another) result(res)
  type(rank2_tensor), intent(in) :: tensor, another
  type(rank2_tensor) :: res
  res%axes = tensor%axes
  res%components = tensor%components + another%components
end function


elemental function add_in_common_axes(tensor, another) result(res)
  !! Returns the sum of two tensors in the cartesian axis system of @p tensor.
  !! 
  !! @p tensor the first tensor defining the axis system.
  !! @p another the second tensor that is rotated before addition.
  !!
  type(rank2_tensor), intent(in) :: tensor, another
  type(rank2_tensor) :: res
  res%axes = tensor%axes
  res = tensor + rotate(another, tensor)
end function


function sum(tensors) result(res)
  type(rank2_tensor), intent(in) :: tensors(:)
  type(rank2_tensor) :: res
  integer :: i
  res%components = 0._dp
  res%axes = tensors(1)%axes
  do i = 2, size(tensors)
    res = res + tensors(i)
  end do
end function

end module
