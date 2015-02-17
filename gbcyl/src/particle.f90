!> Implements the basic functions needed to handle data for simple, 
!! uniaxial particles such as in the Gay-Berne potential.
module particle
use nrtype
use utils
implicit none

!> Holds data of a uniaxial particle, e.g. Gay-Berne 
!! particle. Using the rod variable one can use the same datatype to describe
!! two different particles.
!! 
!! x, y, z are the cartesian positional coordinates of the particle.
!! ux, uy, uz are the components of the unit vector that describes
!! the 
!! orientation of the particle.
!! rod tells if the particle is a rodlike particle or e.g. spherically
!! symmetric.
!!
!! :TODO: Remove rod variable. This module should only describe GB particles.
!!
type particledat
   real(dp) :: x = 0._dp
   real(dp) :: y = 0._dp
   real(dp) :: z = 0._dp
   real(dp) :: ux = 0._dp
   real(dp) :: uy = 0._dp
   real(dp) :: uz = 1._dp
   logical :: rod = .true.
end type particledat

contains

!> Writes @p aparticle to the outputunit @p writeunit.
subroutine writeparticle(writeunit, aparticle)
  integer, intent(in) :: writeunit
  type(particledat), intent(in) :: aparticle
  character(len = 2) :: typeid
  if (aparticle%rod) then
     typeid = 'gb'
  else 
     typeid = 'lj'
  end if
  write(writeunit, '(A,1X,6(' // fmt_char_dp() // ',1X))', advance='no') &
       typeid, aparticle%x, aparticle%y, aparticle%z, aparticle%ux, &
       aparticle%uy, aparticle%uz 
end subroutine writeparticle

!> Reads @p aparticle from @p readunit. If the reading succeeds, 
!! @p ios==0.
subroutine readparticle(readunit, aparticle, ios)
  type(particledat), intent(out) :: aparticle
  integer, intent(in) :: readunit
  integer, intent(out) :: ios
  character(len = 2) :: temp
  read(readunit, fmt=*, iostat=ios) &
       temp, aparticle%x, aparticle%y, aparticle%z, aparticle%ux, &
       aparticle%uy, aparticle%uz 
  if (temp == 'gb') then
     aparticle%rod = .true.
  else if (temp == 'lj') then 
     aparticle%rod = .false. 
  else
     ios = 999
  end if
  if (ios /= 0) backspace readunit
end subroutine readparticle

!> Returns the position of @p aparticle as a vector.
pure function position(aparticle)
  real(dp), dimension(3) :: position
  type(particledat), intent(in) :: aparticle
  position = (/aparticle%x, aparticle%y, aparticle%z/)
end function position

!> Returns the orientation of @p aparticle as a vector. @p vec must
!! be a unit vector.
pure function orientation(aparticle) 
  real(dp), dimension(3) :: orientation
  type(particledat), intent(in) :: aparticle
  orientation = (/aparticle%ux, aparticle%uy, aparticle%uz/)
end function orientation

!> Assigns the position @p vec to @p aparticle.
pure subroutine setposition(aparticle, vec)
  type(particledat), intent(inout) :: aparticle
  real(dp), dimension(3), intent(in) :: vec
  aparticle%x = vec(1)
  aparticle%y = vec(2)
  aparticle%z = vec(3)
end subroutine setposition

!> Assigns the orientation @p vec to @p aparticle.
pure subroutine setorientation(aparticle, vec)
  type(particledat), intent(inout) :: aparticle
  real(dp), dimension(3), intent(in) :: vec
  aparticle%ux = vec(1)
  aparticle%uy = vec(2)
  aparticle%uz = vec(3)
end subroutine setorientation
  
end module
