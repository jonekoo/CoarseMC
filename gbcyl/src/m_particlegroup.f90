module m_particlegroup
use nrtype, only: dp
implicit none

type, abstract :: particlegroup

contains  

!  procedure(pair_reduce), deferred :: reduce_pairs
!  procedure(single_reduce), deferred :: reduce

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

type, extends(particlegroup) :: spheregroup
  type(sphere), allocatable :: spheres(:)
contains
   procedure :: reduce_pairs => reduce_all_pairs
   procedure :: reduce => reduce_no_cutoff
end type spheregroup



type, abstract :: sphere_potential
   contains
     procedure(v_sphere), deferred :: value 
end type sphere_potential

abstract interface
subroutine v_sphere(this, s, res, err)
  use nrtype
  import sphere_potential
  import sphere
  class(sphere_potential), intent(in) :: this
  type(sphere), intent(in), target :: s
  real(dp), intent(out) :: res
  integer, intent(out) :: err
end subroutine v_sphere
end interface



type, abstract :: sphere_pair_potential
   contains
     procedure(v_sphere_pair), deferred :: value 
end type sphere_pair_potential

abstract interface
pure subroutine v_sphere_pair(this, s1, s2, res, err)
  use nrtype
  import sphere_pair_potential
  import sphere
  class(sphere_pair_potential), intent(in) :: this
  type(sphere), intent(in) :: s1, s2
  real(dp), intent(out) :: res
  integer, intent(out) :: err
end subroutine
end interface



contains

!! Self-interaction
pure subroutine reduce_all_pairs(this, potential, res, err)
  class(spheregroup), intent(in) :: this
  class(sphere_pair_potential), intent(in) :: potential
  real(dp), intent(out) :: res
  integer, intent(out) :: err
  real(dp) :: r
  res = 0.0
  err = 0
  !! Tällaiset tarkistukset voisi varmaan tehdä konstruktorissa.
  do i = 2, size(this%spheres)
     do j = 1, i - 1
        call potential%value(this%spheres(i), this%spheres(j), r, err)
        res = res + r
        if (err /= 0) return
     end do
  end do
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





module m_rodgroup
use m_spheregroup

type, extends(sphere) :: rod
  real(dp) :: ux, uy, uz
end type



type, abstract :: rod_potential
   contains
     procedure(v_rod), deferred :: value
end type

abstract interface
subroutine v_rod(this, rod1, res, err)
    use nrtype
    import rod_potential
    import rod
    class(rod_potential), intent(in) :: this
    type(rod), intent(in), target :: rod1
    real(dp), intent(out) :: res
    integer, intent(out) :: err
  end subroutine v_rod
end interface



type, abstract :: rod_pair_potential
  contains
    procedure(v_rod_pair), deferred :: value
end type

abstract interface
  pure subroutine v_rod_pair(this, rod1, rod2, res, err)
    use nrtype
    import rod_pair_potential
    import rod
    class(rod_pair_potential), intent(in) :: this
    type(rod), intent(in) :: rod1, rod2
    real(dp), intent(out) :: res
    integer, intent(out) :: err
  end subroutine v_rod_pair
end interface



type, extends(particlegroup) :: rodgroup
    class(rod_potential), allocatable :: ias(:)
  contains
    procedure :: reduce => reduce_rods
    procedure :: reduce_pairs => reduce_rodpairs
    procedure :: add_interaction => add_rod_ia
end type



contains

subroutine add_rod_ia(this, ia)
  class(rodgroup), intent(inout) :: this
  class(rod_potential), intent(in), target :: ia
  
end subroutine 


pure subroutine reduce_rods(this, potential, res, err)
  class(rodgroup), intent(in) :: this
  class(rod_potential), intent(in) :: potential
  real(dp), intent(out) :: res
  integer, intent(out) :: err
end subroutine

pure subroutine reduce_rodpairs(this, pair_potential, res, err)
  class(rodgroup), intent(in) :: this
  class(rod_pair_potential), intent(in) :: pair_potential
  real(dp), intent(out) :: res
  integer, intent(out) :: err
end subroutine reduce_rodpairs

end module



module m_gblj
use m_rodgroup
use m_spheregroup


type, abstract :: rodsphere_potential
   contains
     procedure(v_rodsphere), deferred :: value
end type

abstract interface
  pure subroutine v_rodsphere(this, ui, rij, energy, overlap)
    use nrtype
    import rodsphere_potential
    class(rodsphere_potential), intent(in) :: this
    real(dp), intent(in) :: ui(3), rij(3)
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
  end subroutine v_rodsphere
end interface

!! This module should be able to return a interaction for rods when the
!! constructor is given spheregroup and vice versa.

type, extends(rod_potential) :: gblj_rods
  type(spheregroup), pointer :: sg => null()
  class(rodsphere_potential), pointer :: potential => null()
  contains
    procedure :: value => value_rods
    final :: finalize_rods
end type

interface gblj_rods
   module procedure construct_rod_ia
end interface


type, extends(sphere_potential) :: gblj_spheres
  type(rodgroup), pointer :: rg => null()
  class(rodsphere_potential), pointer :: potential => null()
  contains
    procedure :: value => value_spheres
    final :: finalize_spheres
end type

interface gblj_spheres
   module procedure construct_sphere_ia
end interface


type, extends(rod_potential) :: gblj_bind_sphere
  type(sphere), pointer :: s => null()
  class(rodsphere_potential), pointer :: potential
contains
  procedure :: value => single_rod_value
end type

type, extends(sphere_potential) :: gblj_bind_rod
  type(rod), pointer :: r => null()
  class(rodsphere_potential), pointer :: potential
contains
  procedure :: value => single_sphere_value
end type


contains 

function construct_rod_ia(potential, sg) result(this)
  class(rodsphere_potential), intent(in), target :: potential
  type(spheregroup), intent(in), target :: sg
  type(gblj_rods) :: this
  this%sg => sg
  this%potential => potential
end function

subroutine value_rods(this, rod1, res, err)
  class(gblj_rods), intent(in) :: this
  type(rod), intent(in), target :: rod1
  real(dp), intent(out) :: res
  integer, intent(out) :: err
  type(gblj_bind_rod) :: pot
  pot%r => rod1 ! bind rod
  pot%potential => this%potential
  call this%sg%reduce(pot, res, err)
end subroutine

pure subroutine single_sphere_value(this, s, res, err)
  class(gblj_bind_rod), intent(in) :: this
  type(sphere), intent(in), target :: s
  real(dp), intent(out) :: res
  integer, intent(out) :: err
  logical :: overlap
  call this%potential%value([this%r%ux, this%r%uy, this%r%uz],&
       [s%x - this%r%x, s%y - this%r%y, s%z - this%r%z], res, overlap)
  if (overlap) then
     res = 1
  else
     res = 0
  end if
end subroutine

subroutine finalize_rods(this)
  type(gblj_rods), intent(inout) :: this
  if (associated(this%sg)) this%sg => null()
end subroutine



function construct_sphere_ia(rg)
  type(rodgroup), intent(in), target :: rg
  type(gblj_spheres) :: construct_sphere_ia
  construct_sphere_ia%rg => rg
end function

subroutine value_spheres(this, s, res, err)
  class(gblj_spheres), intent(in) :: this
  type(sphere), intent(in), target :: s
  real(dp), intent(out) :: res
  integer, intent(out) :: err
  type(gblj_bind_sphere) :: pot
  pot%s => s
  pot%potential => this%potential
  call this%rg%reduce(pot, res, err) 
end subroutine


subroutine finalize_spheres(this)
  type(gblj_spheres), intent(inout) :: this
  if (associated(this%rg)) this%rg => null()
end subroutine


pure subroutine single_rod_value(this, rod1, res, err)
  use gblj
  class(gblj_bind_sphere), intent(in) :: this
  type(rod), intent(in), target :: rod1
  real(dp), intent(out) :: res
  integer, intent(out) :: err
  logical :: overlap
  call this%potential%value([rod1%ux, rod1%uy, rod1%uz], &
       [this%s%x - rod1%x, this%s%y-rod1%y, this%s%z - rod1%z], res, overlap)
  if (overlap) then
     res = 1
  else
     res = 0
  end if
end subroutine


end module m_gblj



