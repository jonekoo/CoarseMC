!> Implements routines for creating trial Monte Carlo translations and
!! rotations for uniaxial particles, such as in the Gay-Berne model.
module particle_mover
use nrtype
use utils
use class_parameterizer
use class_parameter_writer
include 'rng.inc'
implicit none
  
!> The maximum angle for a trial rotation.
real(dp), save :: max_rotation = -1._dp

!> The maximum change in position for a trial translation.
real(dp), save :: max_translation = -1._dp

!> True when this module has been initialized.
logical, save :: is_initialized = .false.

!> Initializes the module. 
interface particlemover_init
  module procedure particlemover_init1, particlemover_init2
end interface

contains

  !> Initializes this module using a parameterizer object.
  !! 
  !! @param reader the parameterizer object which gives (reads) the parameters
  !! to this module.
  !! 
  subroutine particlemover_init1(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'max_translation', max_translation)
    call getparameter(reader, 'max_rotation', max_rotation)
    if(max_translation < 0._dp) stop 'No max_translation given. Stopping.'
    if(max_rotation < 0._dp) then
       !! Set rotation to move end of molecule about as much as translation.
       max_rotation = 2._dp * asin(max_translation / 4.4_dp)
       write(*, *) 'init_particle: Assuming molecule length is 4.4.'
    end if
    is_initialized = .true.
  end subroutine particlemover_init1

  !> Initializes the module by setting the maximum translation and
  !! rotation to @p distance and @p angle, respectively. 
  subroutine particlemover_init2(distance, angle)
    implicit none
    real(dp), intent(in) :: distance 
    real(dp), intent(in) :: angle
    max_rotation = angle
    max_translation = distance
    is_initialized = .true.
  end subroutine

  !> Writes the module parameters. Output unit and format are defined
  !! by @p writer.
  subroutine particlemover_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'Particle moving parameters')
    call writeparameter(writer, 'max_translation', max_translation)
    call writeparameter(writer, 'max_rotation', max_rotation)    
  end subroutine

  !> Creates a new position @p xn, @p yn, @p zn from @p xo, @p yo, @p zo by
  !! applying a random displacement to each of the three coordinates.
  !! @p genstate is the random number generator state. 
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

  !> Creates a new unit vector @p uxn, @p uyn, @p uzn by applying a random
  !! rotation to unit vector @p uxo, @p uyo, @p uzo.
  pure subroutine rotate(uxo, uyo, uzo, uxn, uyn, uzn, genstate)
    real(dp), intent(in) :: uxo, uyo, uzo
    real(dp), intent(out) :: uxn, uyn, uzn
    type(rngstate), intent(inout) :: genstate
    real(dp) :: theta, nx, ny, nz
    real(dp) :: r
    call nvec(nx, ny, nz, genstate)
    call rng(genstate, r)
    theta = (2._dp * r - 1._dp) * max_rotation
    call rotate_vector(uxo, uyo, uzo, nx, ny, nz, theta, uxn, uyn, uzn)
  end subroutine

  !> Sets the maximum translation @p distance and maximum rotation
  !! @p angle for a single particle move.
  subroutine setmaxmoves(distance, angle)
    implicit none
    real(dp), intent(in) :: distance, angle
    max_translation = distance
    max_rotation = angle
  end subroutine

  !> Getter for the maximum translation @p distance and maximum rotation
  !! @p angle for a single particle move.
  subroutine getmaxmoves(distance, angle)
    implicit none
    real(dp), intent(out) :: distance, angle
    distance = max_translation
    angle = max_rotation
  end subroutine 

  !> Returns the maximum possible translation of a particle.
  real(dp) function get_max_translation()
    if (.not. is_initialized) stop 'Trying to access max_translation before module particle is initialized.'
    get_max_translation = max_translation
  end function
    
  !> Generates a random unit vector (@p nx, @p ny, @p nz). @p genstate
  !! is the random number generator state.
  !!
  !! @see Understanding Mol. Sim. 2nd Ed.  Frenkel, Smit p. 578
  !!
  pure subroutine nvec(nx, ny, nz, genstate)
    include 'rng.inc'
    intrinsic sqrt
    double precision, intent(out) :: nx, ny, nz
    type(rngstate), intent(inout) :: genstate
    double precision :: l, u1, u2, s
    double precision :: r
    l = 0.0_dp
    do
       call rng(genstate, r)
       u1 = 1._dp - 2._dp * r
       call rng(genstate, r)
       u2 = 1._dp - 2._dp * r
       l = u1 * u1 + u2 * u2
       if(l <= 1._dp) exit
    end do
    s = 2.0_dp * sqrt(1._dp - l)
    nx = u1 * s
    ny = u2 * s
    nz = 1._dp - 2._dp * l
  end subroutine nvec

end module particle_mover
