program onwallpressure
  use m_particle_factory
  use particle
  use class_poly_box
  use class_parameterizer
  use m_fileunit
  use particlewall
  use num_kind
  use utils
  implicit none
  type(factory) :: coordinatereader
  integer :: coordinateunit
  integer :: ios
  type(particledat), allocatable :: particles(:)
  type(poly_box) :: simbox
  type(parameterizer) :: reader
  character(len=3) :: idchar = '0'
  read(*, *) idchar
  reader = new_parameterizer('inputparameters.'//trim(adjustl(idchar)))
  !! Initialize the modules needed
  call particlewall_init(reader)
  coordinateunit = fileunit_getfreeunit()
  open(unit=coordinateunit, file='configurations.' //trim(adjustl(idchar)), &
       action='READ', status='OLD')
  do
     call factory_readstate(coordinatereader, coordinateunit, simbox, &
          particles, ios)
     if (ios /= 0) then
        exit
     end if
     write(*, *) pressure(simbox, particles)
  end do
  
contains
  
  function pressure(simbox, particles)
    implicit none
    type(poly_box), intent(in) :: simbox
    type(particledat), intent(in) :: particles(:)
    real(dp) :: pressure
    real(dp) :: f(3)
    real(dp) :: radius
    integer :: i
    real(dp) :: sum_f(3)
    real(dp), parameter :: pi = 4._dp * atan(1._dp)
    real(dp) :: r(3)
    radius = getx(simbox) / 2._dp
    pressure = 0._dp
    sum_f = 0._dp
    do i = 1, size(particles)
       f = particlewall_force(particles(i), simbox)
       r = position(particles(i))
       r(3) = 0._dp
       r = r / sqrt(dot_product(r, r))
       pressure = pressure + dot_product(f, r)
    end do
    pressure = pressure / (2 * pi * radius * getz(simbox))
  end function pressure
  
end program onwallpressure
