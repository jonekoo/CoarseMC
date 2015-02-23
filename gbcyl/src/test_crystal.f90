program test_crystal
  use m_crystal
  use particle
  use class_factory
  use class_poly_box
  implicit none
  integer, parameter :: nx = 12, ny = 12, nz = 6
  real(8) :: xs(nx, ny, nz), ys(nx, ny, nz), zs(nx, ny, nz)
  real(8), parameter :: a = 1.0, d = 3.6
  real(8) :: rs(3, nx * ny * nz)
  type(poly_box) :: simbox
  type(particledat) :: particles(nx * ny * nz)
  integer :: i
  integer, parameter :: stdout = 6
  real(8), parameter :: h = sqrt(3.0) / 2 * a
  type(factory) :: f

  call hexagonally_close_packed(nx, ny, nz, a, d, .true., xs, ys, zs)
 
  rs(1, :) = reshape(xs, [nx * ny * nz]) 
  rs(2, :) = reshape(ys, [nx * ny * nz]) 
  rs(3, :) = reshape(zs, [nx * ny * nz]) 

  do i = 1, nx * ny * nz
     call setposition(particles(i), rs(:, i))
  end do

  simbox = new_box(maxval(xs) + a, maxval(ys) + h, maxval(zs) + d)
  call factory_writestate(f, stdout, simbox, particles)
end program

