program create_crystal
  use utils, only: splitstr, join
  use m_crystal
  use particle
  use m_particle_factory
  use class_poly_box
  use iso_fortran_env
  implicit none
  integer :: nx = -1, ny = -1, nz = -1
  real(8), allocatable :: xs(:, :, :), ys(:, :, :), zs(:, :, :) 
  real(8) :: a = 1.0, d = 3.6
  real(8), allocatable :: rs(:, :) !(3, nx * ny * nz)
  type(poly_box) :: simbox
  type(particledat), allocatable :: particles(:)
  integer :: i
  real(8) :: h != sqrt(3.0) / 2 * a
  type(factory) :: f
  integer :: unit = 12
  character(len=80) :: ofile = "configurations.txt"
  character(len=15) :: boxtype = "rectangular"
  logical, allocatable :: mask(:)
  real(dp) :: radius = 9.0,  offset = 0.5

  !!---- Command line parsing ------------------------------------------------
  character(len=1000) :: cmd_line = ""
  namelist /cmd/ nx, ny, nz, a, d, ofile, boxtype, radius, offset
  character(len=*), parameter :: options_str = &
       'ofile="configurations.txt", a=1.0, d=3.6, nx=12, ny=12, nz=6, ' // &
       'boxtype="rectangular", [radius=9.0,] [offset=0.5]'

  call get_command(cmd_line)
  call parse_cmdline(cmd_line)

  if (ofile == "") then
     write(error_unit, *) 'Give ofile!'
     call print_usage
     stop
  else if (nx < 1 .or. ny < 1 .or. nz < 1) then
     write(error_unit, *) 'nx, ny and nz must be > 1!'   
     call print_usage
     stop
  else if (a <= 0.0 .or. d <= 0.0) then
     write(error_unit, *) 'a and d must be > 0!'
     call print_usage
     stop
  end if
  !!--------------------------------------------------------------------------

  h = sqrt(3.0) / 2 * a

  simbox = new_box(boxtype, (nx + 1) * a, (ny + 1) * h, (nz + 1) * d)

  allocate(xs(nx, ny, nz), ys(nx, ny, nz), zs(nx, ny, nz), &
       rs(3, nx * ny * nz), particles(nx * ny * nz))
  call hexagonally_close_packed(nx, ny, nz, a, d, .true., xs, ys, zs)
 
  rs(1, :) = reshape(xs, [nx * ny * nz]) 
  rs(2, :) = reshape(ys, [nx * ny * nz]) 
  rs(3, :) = reshape(zs, [nx * ny * nz]) 

  do i = 1, nx * ny * nz
     call setposition(particles(i), rs(:, i))
  end do

  open(unit=unit, file=ofile, action='write')
  if (gettype(simbox) == 'cylindrical') then
     if (radius > 0) call setx(simbox, radius * 2)
     allocate(mask(size(particles)))
     !call cylinder_mask(particles, getx(simbox) / 2 - offset, mask)
     !mask = .false.
     mask(:) = particles(:)%x**2 + particles(:)%y**2 < (getx(simbox) / 2 - offset)**2
     call factory_writestate(f, unit, simbox, pack(particles, mask))
  else
     call factory_writestate(f, unit, simbox, particles)
  end if
  close(unit)


contains

include 'parse_cmdline.inc'


subroutine print_help
  call print_usage
  write(output_unit, *) 'This program creates a hexagonally close-packed'
  write(output_unit, *) 'crystal and places it in to the simulation box.'
  write(output_unit, *) ''
  write(output_unit, *) 'The options for the program are listed below.'
  write(output_unit, *) ''
  write(output_unit, *) 'ofile=filename'
  write(output_unit, *) '    the output file where the particles and the'
  write(output_unit, *) '     simulation box is written to.'
  write(output_unit, *) ''
  write(output_unit, *) 'a = 1.0'
  write(output_unit, *) '    the distance between nearest neighbours in the'
  write(output_unit, *) '    hexagonal layers in the x,y-plane.'
  write(output_unit, *) ''
  write(output_unit, *) 'boxtype="rectangular"'
  write(output_unit, *) '    options are "cylindrical", "rectangular". You'
  write(output_unit, *) '    may need to escape the quotes with backslash.'
  write(output_unit, *) '    Example: boxtype=\"cylindrical\"'
  write(output_unit, *) ''
  write(output_unit, *) 'd = 3.6'
  write(output_unit, *) '    the distance between layers in the z-direction.'
  write(output_unit, *) ''
  write(output_unit, *) 'nx = 12, ny = 12, nz = 6'
  write(output_unit, *) '    the number of particles in the x-, y- and'
  write(output_unit, *) '    z-directions.'
  write(output_unit, *) ''
  write(output_unit, *) 'radius = 9.0 (optional)'
  write(output_unit, *) '    if boxtype is cylindrical then this option can'
  write(output_unit, *) '    be used to set the radius where the wall of the'
  write(output_unit, *) '    cylinder is.'
  write(output_unit, *) ''
  write(output_unit, *) 'offset = 0.5 (optional)'
  write(output_unit, *) '    if boxtype is cylindrical then this option can'
  write(output_unit, *) '    be used to control how close to the cylinder'
  write(output_unit, *) '    wall the particles are put.'
  write(output_unit, *) ''
end subroutine

subroutine cylinder_mask(theparticles, cutradius, themask)
  type(particledat), intent(in) :: theparticles(:)
  real(dp), intent(in) :: cutradius
  logical, intent(out) :: themask(size(theparticles))
  themask = .false.
  do i = 1, size(theparticles)
     themask(i) = theparticles(i)%x**2+theparticles(i)%y**2 < cutradius**2
  end do
end subroutine

end program

