module m_particlegroup
use nrtype, only: dp
implicit none

type, abstract :: particlegroup

contains  

  procedure(pair_reduce), deferred :: reduce_pairs
  procedure(single_reduce), deferred :: reduce

end type particlegroup


type, abstract :: particle
   real(dp) :: x, y, z
end type particle


abstract interface 
  pure subroutine pair_reduce(this, potential, res, err)
    use nrtype, only: dp
    import particlegroup
    class(particlegroup), intent(in) :: this
    class(*), intent(in) :: potential
    real(dp), intent(out) :: res
    integer, intent(out) :: err
  end subroutine pair_reduce
end interface


abstract interface
  pure subroutine single_reduce(this, potential, res, err)  
    use nrtype, only: dp
    import particlegroup
    class(particlegroup), intent(in) :: this
    class(*), intent(in) :: potential
    real(dp), intent(out) :: res
    integer, intent(out) :: err
   end subroutine single_reduce
end interface

end module m_particlegroup


module m_spheregroup
use nrtype, only: dp
use m_particlegroup


type sphere
  real(dp) :: x, y, z
end type sphere

type :: spherepotential
   contains
     procedure :: value => pot
end type spherepotential

type, extends(particlegroup) :: spheregroup
  type(sphere), allocatable :: spheres(:)
contains
   procedure :: reduce_pairs => reduce_all_pairs
   procedure :: reduce => reduce_no_cutoff
end type spheregroup

contains

pure subroutine pot(this, s1, s2, res, err)
  class(spherepotential), intent(in) :: this
  type(sphere), intent(in) :: s1, s2
  real(dp), intent(out) :: res
  integer, intent(out) :: err
  res = 0.0
  err = 0
end subroutine pot


!! Self-interaction
pure subroutine reduce_all_pairs(this, potential, res, err)
  class(spheregroup), intent(in) :: this
  class(*), intent(in) :: potential
  real(dp), intent(out) :: res
  integer, intent(out) :: err
  real(dp) :: r
  res = 0.0
  err = 0
  !! Tällaiset tarkistukset voisi varmaan tehdä konstruktorissa.
  select type (potential)
  class is (spherepotential)
    do i = 2, size(this%spheres)
       do j = 1, i - 1
          call potential%value(this%spheres(i), this%spheres(j), r, err)
          res = res + r
          if (err /= 0) return
       end do
    end do
  end select
end subroutine


!! Interaction with others
pure subroutine reduce_no_cutoff(this, potential, res, err)
  class(spheregroup), intent(in) :: this
  class(*), intent(in) :: potential
  real(dp), intent(out) :: res
  integer, intent(out) :: err
  res = 0.0
  err = 0
end subroutine


end module m_spheregroup


module m_rodsphere
use m_spheregroup

type, extends(sphere) :: rod
  real(dp) :: ux, uy, uz
end type

type, abstract :: rodsphere_potential   
   type(rod), pointer :: rp => null()
   type(sphere), pointer :: sp => null()
   contains
     procedure(v_rod), deferred :: value_rod
     procedure(v_sphere), deferred :: value_sphere
     generic :: value => value_rod, value_sphere
end type

interface rodsphere_potential
procedure constructor_rod, constructor_sphere
end interface

interface
   subroutine v_rod(this, somerod, res, err)
     use nrtype
     import rodsphere_potential
     import rod
     class(rodsphere_potential), intent(in) :: this
     type(rod), intent(in) :: somerod
     real(dp), intent(out) :: res
     integer, intent(out) :: err
   end subroutine v_rod
end interface

interface
   subroutine v_sphere(this, somesphere, res, err)
     use nrtype
     use m_spheregroup
     import rodsphere_potential
     class(rodsphere_potential), intent(in) :: this
     type(sphere), intent(in) :: somesphere
     real(dp), intent(out) :: res
     integer, intent(out) :: err
   end subroutine v_sphere
end interface



contains

function constructor_rod(somerod)
  type(rod), intent(in), target :: somerod
end function constructor_rod

function constructor_sphere(somesphere)
  type(sphere), intent(in) :: somesphere
end function constructor_sphere



end module m_rodsphere

