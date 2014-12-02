!! :TODO: Think what functions and subroutines could be made elemental.
module particle_mover
  use particle
  use nrtype
  use utils
  use class_parameterizer
  use class_parameter_writer
  include 'rng.inc'
  implicit none
  private

  public :: initparticle
  public :: particle_writeparameters
  public :: getmaxmoves
  public :: move
  public :: setmaxmoves
  public :: get_max_translation
  
  real(dp), save :: dthetamax = -1._dp
  real(dp), save :: max_translation = -1._dp
  logical, save :: is_initialized = .false.

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
    if(dthetamax < 0._dp) then
      !! Set rotation to move end of molecule about as much as translation.
      dthetamax = 2._dp * asin(max_translation / 4.4_dp)
      write(*, *) 'init_particle: Assuming molecule length is 4.4.'
    end if
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
  
pure subroutine nvec(nx, ny, nz, genstate)
  !!
  !! Generates a random unit vector (nx, ny, nz). 
  !!
  !! @see Understanding Mol. Sim. 2nd Ed.  Frenkel, Smit p. 578
  !!
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
end subroutine

end module particle_mover
