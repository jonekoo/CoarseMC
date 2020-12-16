!> Contains procedures for creating positions of molecules on a periodic
!! crystal lattice.
module m_crystal
implicit none

contains

!> Returns positions @p xs, @p ys, @p zs for a hexagonally close-packed
!! crystal with @p nx, @p ny and @p nz particles in x-, y- and
!! z-directions, respectively. The distance between particles in the
!! hexagonal layers is @p a and the distance between layers is @p d.
!! If @p center is .true. the coordinates are transformed so that the
!! origin is in the middle of the crystal. Otherwise all coordinates are
!! non-negative.
!! 
!! @post xs(i), ys(i), zs(i) is the position of the i:th point (or
!! particle) in the crystal.
!!
subroutine hexagonally_close_packed(nx, ny, nz, a, d, center, xs, ys, zs)
  integer, intent(in) :: nx, ny, nz
  real(8), intent(in) :: a, d
  logical, intent(in) :: center
  real(8), intent(out) :: xs(nx, ny, nz), ys(nx, ny, nz), zs(nx, ny, nz)
  integer :: i
  real(8) :: h !< The distance between adjacent rows in the same layer
  real(8) :: m !< The distance between corresponding rows in adjacent
               !< layers.
  xs = 0.0
  ys = 0.0
  zs = 0.0
  !! First create coordinates, then center.
  !! Create first row in first layer
  xs(:, 1, 1) = [(i * a, i = 0, nx - 1)]
  !! Create other rows in the first layer
  do i = 3, ny, 2
     xs(:, i, 1) = xs(:, 1, 1)
  end do
  do i = 2, ny, 2
     xs(:, i, 1) = xs(:, 1, 1) + 0.5 * a
  end do

  h = sqrt(3.0) / 2 * a
  do i = 1, ny
     ys(:, i, 1) = (i - 1) * h
  end do

  !! Add distances between layers
  do i = 1, nz
     zs(:, :, i) = (i - 1) * d
  end do

  !! For every other layer, x and y are the same.
  !! Odd layers:
  do i = 1, nz, 2
     xs(:, :, i) = xs(:, :, 1)
     ys(:, :, i) = ys(:, :, 1)
  end do
  m = a / (2 * sqrt(3.0))
  !! Even layers:
  do i = 2, nz, 2
     xs(:, :, i) = xs(:, :, 1) + 0.5 * a
     ys(:, :, i) = ys(:, :, 1) + m
  end do
  if (center) then
     xs = xs - maxval(xs) / 2.0
     ys = ys - maxval(ys) / 2.0
     zs = zs - maxval(zs) / 2.0
  end if
end subroutine hexagonally_close_packed

end module m_crystal


