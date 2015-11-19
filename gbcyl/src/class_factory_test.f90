!> Unit tests for the module class_factory.
module class_factory_test
use ftnunit
use class_poly_box
use particle
use class_factory
use num_kind
implicit none

contains

subroutine test_all
  call test(test_writeread, "Writing a configuration to a file and reading it back.") 
end subroutine

subroutine test_writeread
real(dp) :: margin = 1.e-9_dp
integer, parameter :: unit = 999
integer :: ios
type(poly_box) :: boxwritten
type(poly_box) :: boxwritten2
type(particledat), dimension(100) :: particleswritten
type(particledat), dimension(100) :: particleswritten2
type(poly_box) :: boxread
type(particledat), allocatable :: particlesread(:)
real(dp), dimension(100) :: numbers
real(dp), parameter :: diameter = 11.11_dp
real(dp), parameter :: height = 22.22_dp
type(factory) :: thefactory

call random_number(numbers)
particleswritten%x = numbers
call random_number(numbers)
particleswritten%y = numbers
call random_number(numbers)
particleswritten%z = numbers

boxwritten = new_cylinder(diameter, height)

open(unit, status='SCRATCH')
call factory_writestate(thefactory, unit, boxwritten, particleswritten)

boxwritten2 = new_box(height, height, diameter)

!! Write another configuration to file
call random_number(numbers)
particleswritten2%x = numbers
call random_number(numbers)
particleswritten2%y = numbers
call random_number(numbers)
particleswritten2%z = numbers
call factory_writestate(thefactory, unit, boxwritten2, particleswritten2)

!! Read and check the first configuration
rewind unit
call factory_readstate(thefactory, unit, boxread, particlesread, ios) 
!! Check for errors in reading.
call assert_equal(0, ios, "Box could not be read")
!! Check box type
call assert_true('cylindrical' == gettype(boxread), "Box read is not a cylinder")
!! Check box dimensions.
call assert_comparable(diameter, getx(boxread), margin, "Dimensions in x direction &
not comparable")
call assert_comparable(diameter, gety(boxread), margin, "Dimensions in y direction &
not comparable")
call assert_comparable(height, getz(boxread), margin, "Dimensions in z direction &
not comparable")
!! Check particle coordinates.
call assert_comparable(particleswritten%x, particlesread%x, margin, "Particles' &
x-coordinates not comparable.")
call assert_comparable(particleswritten%y, particlesread%y, margin, "Particles' &
y-coordinates not comparable.")
call assert_comparable(particleswritten%z, particlesread%z, margin, "Particles' &
z-coordinates not comparable.")

!! Read and check the second configuration
call factory_readstate(thefactory, unit, boxread, particlesread, ios)
call assert_equal(0, ios, "Box could not be read")
!! Check box dimensions.
call assert_comparable(height, getx(boxread), margin, "Dimensions in x direction &
not comparable")
call assert_comparable(height, gety(boxread), margin, "Dimensions in y direction &
not comparable")
call assert_comparable(diameter, getz(boxread), margin, "Dimensions in z direction &
not comparable")
!! Check particle coordinates.
call assert_comparable(particleswritten2%x, particlesread%x, margin, "Particles' &
x-coordinates not comparable.")
call assert_comparable(particleswritten2%y, particlesread%y, margin, "Particles' &
y-coordinates not comparable.")
call assert_comparable(particleswritten2%z, particlesread%z, margin, "Particles' &
z-coordinates not comparable.")
close(unit)

end subroutine

end module
