!! :TODO: Think what functions and subroutines could be made elemental.
module particle
  use nrtype
  use utils
  implicit none
  private

  public :: particledat
  !public :: createparticle
  public :: position
  public :: orientation
  public :: setposition
  public :: setorientation
  public :: read
  public :: write
  
  !> Holds data of a uniaxial particle, e.g. Gay-Berne 
  !! particle. Using the rod variable one can use the same datatype to describe
  !! two diffrent particles.
  !! 
  !! x, y, z are the cartesian positional coordinates of the particle.
  !! ux, uy, uz are the components of the vector that describes the 
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
  end type

  interface read
    module procedure readparticle
  end interface
 
  interface write
    module procedure writeparticle
  end interface

  contains

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
  end subroutine

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
  end subroutine

  pure function position(aparticle)
    real(dp), dimension(3) :: position
    type(particledat), intent(in) :: aparticle
    position = (/aparticle%x, aparticle%y, aparticle%z/)
  end function
  
  pure function orientation(aparticle) 
    real(dp), dimension(3) :: orientation
    type(particledat), intent(in) :: aparticle
    orientation = (/aparticle%ux, aparticle%uy, aparticle%uz/)
  end function 

  pure subroutine setposition(aparticle, vec)
    type(particledat), intent(inout) :: aparticle
    real(dp), dimension(3), intent(in) :: vec
    aparticle%x = vec(1)
    aparticle%y = vec(2)
    aparticle%z = vec(3)
  end subroutine

  pure subroutine setorientation(aparticle, vec)
    type(particledat), intent(inout) :: aparticle
    real(dp), dimension(3), intent(in) :: vec
    aparticle%ux = vec(1)
    aparticle%uy = vec(2)
    aparticle%uz = vec(3)
  end subroutine

end module
