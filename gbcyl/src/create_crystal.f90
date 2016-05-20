program create_crystal
  use utils, only: splitstr, join
  use m_crystal
  use m_point, only: point, particlearray_to_json
  use class_poly_box
  use iso_fortran_env
  use json_module
  use num_kind, only: dp
  use m_rod, only: rod
  use m_doublerod, only: doublerod
  implicit none
  integer :: nx = -1, ny = -1, nz = -1
  real(8), allocatable :: xs(:, :, :), ys(:, :, :), zs(:, :, :) 
  real(8) :: a = 1.0, d = 3.6
  real(8), allocatable :: rs(:, :) !(3, nx * ny * nz)
  type(poly_box) :: simbox
  class(point), allocatable, target :: particles(:), temp(:)
  integer :: i, j
  real(8) :: h != sqrt(3.0) / 2 * a
  character(len=80) :: ofile = "geometry.json"
  character(len=15) :: boxtype = "rectangular"
  character(len=15) :: particletype = "rod"
  character(len=20) :: groupname = "gb"
  logical, allocatable :: mask(:)
  real(dp) :: radius = 9.0,  offset = 0.5

  !!---- Command line parsing ------------------------------------------------
  character(len=1000) :: cmd_line = ""
  namelist /cmd/ nx, ny, nz, a, d, ofile, boxtype, radius, offset, particletype, &
       groupname
  character(len=*), parameter :: options_str = &
       'ofile="geometry.json", a=1.0, d=3.6, nx=12, ny=12, nz=6, ' // &
       'boxtype="rectangular", particletype="rod", groupname="gb", ' // &
       '[radius=9.0,] [offset=0.5]'
  type(json_value), pointer :: json_val => null(), box_json => null(), &
       groups_json => null(), groups_json_element => null()
  integer, allocatable :: indices(:)

  call get_command(cmd_line)
  call parse_cmdline(cmd_line)

  if (ofile == "") then
     write(error_unit, *) 'Give ofile!'
     call print_usage
     stop
  else if (nx < 1 .or. ny < 1 .or. nz < 1) then
     write(error_unit, *) 'nx, ny and nz must be > 0!'   
     call print_usage
     stop
  else if (a <= 0.0 .or. d <= 0.0) then
     write(error_unit, *) 'a and d must be > 0!'
     call print_usage
     stop
  else if (groupname == "") then
     write(error_unit, *) 'groupname must not be empty!'
     call print_usage
     stop
  end if
  !!--------------------------------------------------------------------------

  h = sqrt(3.0) / 2 * a

  simbox = new_box(boxtype, (nx + 1) * a, (ny + 1) * h, (nz + 1) * d)

  allocate(xs(nx, ny, nz), ys(nx, ny, nz), zs(nx, ny, nz), &
       rs(3, nx * ny * nz))

  call hexagonally_close_packed(nx, ny, nz, a, d, .true., xs, ys, zs)
 
  rs(1, :) = reshape(xs, [nx * ny * nz]) 
  rs(2, :) = reshape(ys, [nx * ny * nz]) 
  rs(3, :) = reshape(zs, [nx * ny * nz]) 

  select case (particletype)
  case ("rod")
    allocate(particles(nx * ny * nz), source=rod())
  case ("point")
    allocate(particles(nx * ny * nz), source=point())
  case ("doublerod")
    allocate(particles(nx * ny * nz), source=doublerod())
  case default
    write(error_unit, *) 'ERROR: Unknown particle type ', particletype
    stop 'create_crystal is unable to continue.'
  end select

  do i = 1, nx * ny * nz
     call particles(i)%set_position(rs(:, i))
  end do

  call json_initialize()

  call json_create_object(json_val, 'crystal')
  call json_create_array(groups_json, 'particle_groups')

  call json_create_object(groups_json_element, '')
  if (gettype(simbox) == 'cylindrical') then
     if (radius > 0) call setx(simbox, radius * 2)
     allocate(mask(size(particles)))
     mask(:) = particles(:)%x**2 + particles(:)%y**2 < &
          (getx(simbox) / 2 - offset)**2

     allocate(temp(count(mask)), source=particles(1:count(mask)))
     j = 0
     do i = 1, size(mask)
        if (mask(i)) then
           j = j + 1
           call temp(j)%downcast_assign(particles(i))
        end if
     end do
     
     call particlearray_to_json(groups_json_element, temp)
  else
     call particlearray_to_json(groups_json_element, particles)
  end  if
  call json_add(groups_json_element, 'name', groupname)
  call json_add(groups_json, groups_json_element)

  call json_create_object(box_json, 'box')
  call simbox%to_json(box_json)
  call json_add(json_val, box_json)

  call json_add(json_val, groups_json)
  call json_print(json_val, ofile)

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
  write(output_unit, *) 'particletype="rod"'
  write(output_unit, *) '    options are "rod", "point". You may need to'
  write(output_unit, *) '    escape the quotes with backslash. Example:'
  write(output_unit, *) '    particletype=\"point\"'
  write(output_unit, *) ''
end subroutine

subroutine cylinder_mask(theparticles, cutradius, themask)
  type(point), intent(in) :: theparticles(:)
  real(dp), intent(in) :: cutradius
  logical, intent(out) :: themask(size(theparticles))
  themask = .false.
  do i = 1, size(theparticles)
     themask(i) = theparticles(i)%x**2+theparticles(i)%y**2 < cutradius**2
  end do
end subroutine

end program

