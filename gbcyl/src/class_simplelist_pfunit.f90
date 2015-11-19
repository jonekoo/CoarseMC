module class_simplelist_pfunit
use pfunit
use particle
use class_poly_box
use class_simplelist
use num_kind
implicit none

contains

subroutine test_new_simplelist
  integer, parameter :: n=3
  type(particledat), dimension(n) :: particles
  type(poly_box) :: simbox
  type(simplelist) :: sl
  integer :: i
  real(dp), parameter :: minlength=10._dp
  simbox = new_cylinder(2._dp*(minlength+tiny(minlength)), 2._dp*(minlength+tiny(minlength)))
  do i=1,n
    particles(i)%x = i*(-minlength*0.25_dp)
    particles(i)%y = i*(-minlength*0.25_dp)
    particles(i)%z = i*(-minlength*0.25_dp)
  end do
  !! Create new list 
  call new_simplelist(sl, simbox, particles, minlength, min_boundary_width=0._dp)
  !! Tests for lists contents
  call AssertEqual(n, sl%counts(0,0,0), "Number of particle indices in the&
  & cell 1,1,1 does not match the total particle number.")
  call AssertEqual((/1,2,3/), sl%indices(1:3,0,0,0), "Some particles did not&
  & end up to the right cell.")
  call AssertEqual(n, sum(sl%counts), "Total count of particle indices in cell&
  & list does not match particle count.")
  call simplelist_deallocate(sl)

  !! Create new list with more cells. 
  do i=1,n
    particles(i)%x = i*(minlength-1e-9)*0.25_dp
    particles(i)%y = i*(minlength-1e-9)*0.25_dp
    particles(i)%z = i*(minlength-1e-9)*0.25_dp
  end do
  call new_simplelist(sl, simbox, particles, minlength, min_boundary_width=0._dp)
  call AssertEqual(2, sl%counts(2,2,2), "More or less than two particles ended up in&
  & the cell 2,2,2")
  call AssertEqual(1, sl%counts(3,3,3), "More or less than one particle ended up in&
  & the cell 3,3,3")
  call AssertEqual((/1,2/), sl%indices(1:2,2,2,2), "Particles 1,2 did not end&
  & up as first particles in the cell 2,2,2")
  call AssertEqual((/3,0/), sl%indices(1:2,3,3,3), "Particle 3 did not end up as&
  & the first particle in cell 3,3,3")
  call AssertEqual(n, sum(sl%counts), "Total particle count does&
  & not match the total count of particle indices in the cell list.")
  call simplelist_deallocate(sl)

  !! Cases where the creation might break:
  !! 1. Only one cell in some dimension
  call new_simplelist(sl, simbox, particles, minlength, min_boundary_width=0._dp)
  call AssertEqual(n, sl%counts(0,0,0), "Number of particle indices in the&
  & cell 1,1,1 does not match the total particle number.")
  call AssertEqual((/1,2,3/), sl%indices(1:3,0,0,0), "Some particles did not&
  & end up to the right cell.")
  call AssertEqual(n, sum(sl%counts), "Total count of particle indices in cell&
  & list does not match particle count.")
  call simplelist_deallocate(sl)
  !! 2. iseven=.true. but only one cell can fit with the given minlength
  !! What should happen? Failure? Warning? Silence and only one cell created?
  !! Failure. But how is this tested? There could be a return code but the 
  !! user does not have to handle a return code, so this is not a safe solution. 
  !! Let's just make it stop.
  !!call simplelist_init(minlength=11._dp, iseven=.true.)
  !!call new_simplelist(sl, simbox, particles) !! Will cause a Fortran STOP
end subroutine

subroutine test_updateall
  integer, parameter :: n=2
  type(particledat), dimension(n) :: particles
  type(poly_box) :: simbox
  type(simplelist) :: sl
  real(dp), parameter :: minlength=5._dp
  real(dp), parameter :: boxside=4._dp*(minlength+tiny(minlength))
  simbox=new_box(boxside,boxside,boxside)
  !! Put one particle in (1,1,1) and another in (3,3,3)
  particles(1)%x=-1.5_dp*minlength
  particles(1)%y=-1.5_dp*minlength
  particles(1)%z=-1.5_dp*minlength
  particles(2)%x=0.5_dp*minlength
  particles(2)%y=0.5_dp*minlength
  particles(2)%z=0.5_dp*minlength
  !! Make a list of eight cells
  call new_simplelist(sl, simbox, particles, minlength, min_boundary_width=0._dp)
  call AssertEqual(n,sum(sl%counts), "Count of particle indices in cell list&
  & does not match particle count.")
  call AssertEqual(1,sl%counts(0,0,0), "More or less than one particle in cell&
  & 1,1,1")
  call AssertEqual(1,sl%counts(2,2,2), "More or less than one particle in cell&
  & 2,2,2")
  !! Move first particle to (0,1,0)
  particles(1)%y=-0.5_dp*minlength
  !! Move second particle to (2,2,1)
  particles(2)%z=-0.5_dp*minlength
  call simplelist_update(sl, simbox, particles)
  !! Check positions.
  call AssertEqual(n,sum(sl%counts), "Count of particle indices in cell list&
  & does not match particle count after update.")
  call AssertEqual(1,sl%counts(0,1,0), "More or less than one particle in cell&
  & 0,1,0")
  call AssertEqual(1,sl%counts(2,2,1), "More or less than one particle in cell&
  & 2,2,1")
  call simplelist_deallocate(sl)
end subroutine

subroutine test_nbrmask
  integer, parameter :: n=3
  real(dp), parameter :: minlength=5._dp
  real(dp), parameter :: boxside=nearest(4._dp*minlength, 1._dp)
  type(poly_box) :: simbox
  type(particledat), dimension(n) :: particles
  type(simplelist) :: sl
  logical, dimension(n) :: mask
  !! set up
  simbox=new_cylinder(boxside,boxside)
  !! Test periodicity
  particles(1)%x=-1.5_dp*minlength
  particles(1)%y=-1.5_dp*minlength
  particles(1)%z=-1.5_dp*minlength
  !! Particle 2 is a periodic neighbour of particle 1
  particles(2)%x=-0.5_dp*minlength
  particles(2)%y=-0.5_dp*minlength
  particles(2)%z= 1.5_dp*minlength
  !! Particle 3 is an ordinary neighbour of particle 2 but not a neighbour of
  !! particle 1
  particles(3)%x=-1.5_dp*minlength
  particles(3)%y=-0.5_dp*minlength
  particles(3)%z=0.5*minlength
  call new_simplelist(sl, simbox,particles, minlength, min_boundary_width=0._dp)
  call AssertEqual((/4,4,4/), shape(sl%counts), "The list has wrong dimensions.")
  call simplelist_nbrmask(sl,simbox, 1, mask)
  call AssertFalse(mask(1), "Periodicity test failed for particle 1")
  call AssertTrue(mask(2), "Periodicity test failed for particle 1")
  call AssertFalse(mask(3), "Periodicity test failed for particle 1")
  call simplelist_nbrmask(sl,simbox, 2, mask)
  call AssertTrue(mask(1), "Periodicity test failed for particle 2")
  call AssertFalse(mask(2), "Periodicity test failed for particle 2")
  call AssertTrue(mask(3), "Periodicity test failed for particle 2")
  call simplelist_nbrmask(sl,simbox, 3, mask)
  call AssertFalse(mask(1), "Periodicity test failed for particle 3")
  call AssertTrue(mask(2), "Periodicity test failed for particle 3")
  call AssertFalse(mask(3), "Periodicity test failed for particle 3")
  
  call simplelist_deallocate(sl)
end subroutine

end module
