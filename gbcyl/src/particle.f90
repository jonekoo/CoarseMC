!Partikkelidatan ja siihen kohdistuvien operaatioiden
!määrittelyt

!! :TODO: Think what functions and subroutines could be made elemental.

module particle
  use nrtype
  use utils
  use class_parameterizer
  use class_parameter_writer
  include 'rng.inc'
  implicit none
  private

  public :: particledat
  public :: initparticle
  public :: createparticle
!  public :: binwrite
  public :: write
  public :: read
  public :: position
  public :: orientation
  public :: particle_writeparameters
  public :: getmaxmoves
  public :: move
  public :: setmaxmoves
  public :: setposition
  public :: setorientation
  public :: get_max_translation
  
  !> Holds data of a cylindrically symmetric particle, e.g. Gay-Berne 
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

  real(dp), save :: dthetamax = -1._dp
  real(dp), save :: max_translation = -1._dp
  logical, save :: is_initialized = .false.
  !!character(len = *), parameter :: typeid = 'gb' 

  interface read
    module procedure readparticle
  end interface
 
  interface write
    module procedure writeparticle
  end interface

!  interface binwrite
!    module procedure binwriteparticle
!  end interface

  !! Generic interface for initialization of the module..
  !! 
  interface initparticle
    module procedure init_particleold, init_particleparameterizer
  end interface

  !! Interface for creating trial moves for the particle.. 
  !!
  interface move
    module procedure move1, move2
  end interface

  contains

  !! Initializes this module using a parameterizer object.
  !! 
  !! @p reader the parameterizer object which gives (reads) the parameters to
  !! this module.
  !! 
  subroutine init_particleparameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'max_translation', max_translation)
    call getparameter(reader, 'max_rotation', dthetamax)
    if(max_translation < 0._dp) stop 'No max_translation given. Stopping.'
    if(dthetamax < 0._dp) stop 'No max_rotation given. Stopping.'
    is_initialized = .true.
  end subroutine

  !! Initializes the module parameters max_translation and dthetamax.
  !!
  subroutine init_particleold(maxTranslation, maxRotation)
    implicit none
    real(dp), intent(in) :: maxTranslation 
    real(dp), intent(in) :: maxRotation
    dthetamax = maxRotation
    max_translation = maxTranslation
    is_initialized = .true.
  end subroutine

  !! Saves the module parameters. The method for saving is defined by @p writer.
  !! 
  !! @p writer the object that implements the saving method for parameters.
  !!
  subroutine particle_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'Particle parameters')
    call writeparameter(writer, 'max_translation', max_translation)
    call writeparameter(writer, 'max_rotation', dthetamax)    
  end subroutine

  !! Currently not in use!!!
  function createparticle(particlestring) result(pp)
    type(particledat) :: pp
    character(len = *), intent(in) :: particlestring
    character(len = 3) :: typeid
    read(particlestring, *) typeid
    if ('gb' == typeid) then
      read(particlestring, *) typeid, pp%x, pp%y, pp%z, pp%ux, pp%uy, pp%uz
      pp%rod = .true.
    else if ('lj' == typeid) then
      read(particlestring, *) typeid, pp%x, pp%y, pp%z
      pp%rod = .false.
    else
      stop 'Error: Could not read particle from string. Stopping.'
    end if
  end function

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

  !subroutine binwriteparticle(writeunit, aparticle)
  !  integer, intent(in) :: writeunit
  !  type(particledat), intent(in) :: aparticle
  !  write(writeunit) &
  !  typeid, aparticle%x, aparticle%y, aparticle%z, aparticle%ux, &
  !  aparticle%uy, aparticle%uz 
  !end subroutine

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

  pure subroutine move1(aparticle, genstate)    
    type(particledat), intent(inout) :: aparticle
    type(rngstate), intent(inout) :: genstate
    type(particledat) :: temp
    call move(aparticle, temp, genstate)
    aparticle = temp
  end subroutine

  pure subroutine move2(oldp, newp, genstate)
    type(particledat), intent(in) :: oldp
    type(particledat), intent(out) :: newp
    type(rngstate), intent(inout) :: genstate
    newp = oldp
    call transmove(oldp%x,oldp%y,oldp%z,newp%x,newp%y,newp%z, genstate)
    if (oldp%rod) then
      call rotate(oldp%ux,oldp%uy,oldp%uz,newp%ux,newp%uy,newp%uz, genstate)
    end if
  end subroutine
  
  pure subroutine transmove(xo, yo, zo, xn, yn, zn, genstate)
    real(dp), intent(in) :: xo, yo, zo
    real(dp), intent(out) :: xn, yn, zn
    type(rngstate), intent(inout) :: genstate
    real(dp) :: max1d 
    real(dp) :: r
    max1d = max_translation/sqrt(3._dp)
    call rng(genstate, r)
    xn = xo + (2._dp * r - 1._dp) * max1d
    call rng(genstate, r)
    yn = yo + (2._dp * r - 1._dp) * max1d
    call rng(genstate, r)
    zn = zo + (2._dp * r - 1._dp) * max1d
  end subroutine

  pure subroutine rotate(uxo, uyo, uzo, uxn, uyn, uzn, genstate)
    real(dp), intent(in) :: uxo,uyo,uzo
    real(dp), intent(out) :: uxn,uyn,uzn
    type(rngstate), intent(inout) :: genstate
    real(dp) :: theta, nx, ny, nz
    real(dp) :: r
    call nvec(nx, ny, nz, genstate)
    call rng(genstate, r)
    theta = (2._dp * r - 1._dp) * dthetamax
    call XVEC2(uxo, uyo, uzo, nx, ny, nz, theta, uxn, uyn, uzn)
  end subroutine

  subroutine setmaxmoves(distance, angle)
    implicit none
    real(dp), intent(in) :: distance, angle
    max_translation = distance
    dthetamax = angle
  end subroutine

  subroutine getmaxmoves(distance, angle)
    implicit none
    real(dp), intent(out) :: distance, angle
    distance = max_translation
    angle = dthetamax
  end subroutine 

  !> Returns the maximum translation of a particle. Can be used by other modules
  !! to get the private max_translation value.
  !!
  real(dp) function get_max_translation()
    if (.not. is_initialized) stop 'Trying to access max_translation before module particle is initialized.'
    get_max_translation = max_translation
  end function
  
end module
