module m_particle
  use iso_fortran_env
  
  type, abstract :: particle
   contains
     procedure(particle_energy), deferred :: energy 
  end type particle

  interface
     subroutine particle_energy(this, res, err)
       import particle, REAL64
       class(particle), intent(inout) :: this
       real(REAL64), intent(out) :: res
       integer, intent(out) :: err
     end subroutine particle_energy
  end interface
  
end module m_particle


module m_sphere
  use iso_fortran_env
  use m_particle
  
  type, abstract, extends(particle) :: sphere
     real(REAL64) :: position(3)
  end type sphere
  
#define MIXINPARAM sphere
#include "m_particlemixin-def.inc"
contains
#include "m_particlemixin-proc.inc"
#undef MIXINPARAM
 
end module m_sphere



module m_rod
  use iso_fortran_env
  use m_particle
  implicit none

  type, abstract, extends(particle) :: rod
     real(REAL64) :: orientation(3) = [0, 0, 1]
  end type rod

#define MIXINPARAM rod
#include "m_particlemixin-def.inc"
contains
#include "m_particlemixin-proc.inc"
#undef MIXINPARAM

end module m_rod

module m_rodext
  use iso_fortran_env
  use m_sphere
  use m_rod, rodmixin => particlemixin
  implicit none
  
  type, extends(rodmixin) :: rodext
   contains
     procedure :: reduce_with_sphere
  end type rodext

  type, abstract :: rodsphere_potential
   contains
     procedure(rodsphere_value), deferred :: value
  end type rodsphere_potential

  interface
     subroutine rodsphere_value(this, the_rod, the_sphere, res, err)
       import rodsphere_potential, rod, sphere, REAL64
       class(rodsphere_potential), intent(in) :: this
       class(rod), intent(in) :: the_rod
       class(sphere), intent(in) :: the_sphere
       real(REAL64), intent(out) :: res
       integer, intent(out) :: err
     end subroutine rodsphere_value
  end interface
  
contains

  subroutine reduce_with_sphere(this, potential, the_sphere, res, err)
    class(rodext), intent(in) :: this
    class(rodsphere_potential), intent(in) :: potential
    class(sphere), intent(in) :: the_sphere
    real(REAL64), intent(out) :: res
    integer, intent(out) :: err
    call potential%value(the_rod=this, the_sphere=the_sphere, res=res, err= err)
  end subroutine reduce_with_sphere
  
end module m_rodext

!#define MIXINPARAM rod
!module m_rodmixin
!  use m_rod
!  use iso_fortran_env
!#include "m_particlemixin-def.inc"
!contains
!#include "m_particlemixin-proc.inc"
!end module m_rodmixin
!#undef MIXINPARAM

!#define MIXINPARAM sphere
!module m_spheremixin
!  use m_sphere
!  use iso_fortran_env
!#include "m_particlemixin-def.inc"
!contains
!#include "m_particlemixin-proc.inc"
!end module m_spheremixin
!#undef MIXINPARAM

program test_mixin
  use m_rod, rodmixin => particlemixin
  use m_sphere, spheremixin => particlemixin

  type(rodmixin) :: rm
  type(spheremixin) :: sm
  real(REAL64) :: res
  integer :: err
  rm%name = 'rod'
  call rm%energy(res, err)
  sm%name = 'sphere'
  call sm%energy(res, err)
  
end program test_mixin
