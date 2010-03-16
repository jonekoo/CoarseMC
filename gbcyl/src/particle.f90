!Partikkelidatan ja siihen kohdistuvien operaatioiden
!m‰‰rittelyt

!! :TODO: Think what functions and subroutines could be made elemental.

module particle
  use nrtype
  use mtmod
  use utils
  use class_parameterizer
  use class_parameter_writer
  implicit none
  private

  public :: particledat
  public :: initparticle
  public :: new_particle
  public :: create_particle
  public :: write_particle
  public :: position
  public :: orientation
  public :: particle_write_parameters
  public :: getmaxmoves
  public :: move
  public :: setmaxmoves
  public :: set_position
  public :: max_trans
  
  type particledat
     real(dp) :: x, y, z, ux, uy, uz
     logical :: rod
  end type particledat

  real(dp), save :: dthetamax = -1._dp
  real(dp), save :: maxdr = -1._dp
  character(len = *), parameter :: type_id_ = 'gb'

  interface initparticle
    module procedure init_particle_old, init_particle_parameterizer
  end interface

  interface max_trans
    module procedure max_trans_f
  end interface

  interface move
    module procedure move1, move2
  end interface

  contains

  subroutine init_particle_parameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call get_parameter(reader, 'max_translation', maxdr)
    call get_parameter(reader, 'max_rotation', dthetamax)
    if(maxdr < 0._dp) stop 'No max_translation given. Stopping.'
    if(dthetamax < 0._dp) stop 'No max_rotation given. Stopping.'
  end subroutine

  subroutine init_particle_old(maxTranslation, maxRotation)
    implicit none
    real(dp), intent(in) :: maxTranslation 
    real(dp), intent(in) :: maxRotation
    dthetamax = maxRotation
    maxdr = maxTranslation
  end subroutine

  subroutine particle_write_parameters(writer)
    type(parameter_writer), intent(in) :: writer
    call write_comment(writer, 'Particle parameters')
    call write_parameter(writer, 'max_translation', maxdr)
    call write_parameter(writer, 'max_rotation', dthetamax)    
  end subroutine

  function new_particle()
    type(particledat) :: new_particle
    new_particle%x = 0._dp
    new_particle%y = 0._dp
    new_particle%z = 0._dp
    new_particle%ux = 0._dp
    new_particle%uy = 0._dp
    new_particle%uz = 1._dp
    new_particle%rod = .true.
  end function

  function create_particle(particle_string) result(pp)
    type(particledat) :: pp
    character(len = *), intent(in) :: particle_string
    character(len = 3) :: type_id
    read(particle_string, *) type_id
    if ('gb' == type_id) then
      read(particle_string, *) type_id, pp%x, pp%y, pp%z, pp%ux, pp%uy, pp%uz
      pp%rod = .true.
    else
      stop 'Error: Could not read particle from string. Stopping.'
    end if
  end function

  subroutine write_particle(write_unit, a_particle)
    integer, intent(in) :: write_unit
    type(particledat), intent(in) :: a_particle
    write(write_unit, '(A, 6' // fmt_char_dp() // ')', advance='no') &
    type_id_, a_particle%x, a_particle%y, a_particle%z, a_particle%ux, &
    a_particle%uy, a_particle%uz 
  end subroutine

  pure function position(a_particle)
    real(dp), dimension(3) :: position
    type(particledat), intent(in) :: a_particle
    position = (/a_particle%x, a_particle%y, a_particle%z/)
  end function
  
  pure function orientation(a_particle) 
    real(dp), dimension(3) :: orientation
    type(particledat), intent(in) :: a_particle
    orientation = (/a_particle%ux, a_particle%uy, a_particle%uz/)
  end function 

  pure subroutine set_position(a_particle, vec)
    type(particledat), intent(inout) :: a_particle
    real(dp), dimension(3), intent(in) :: vec
    a_particle%x = vec(1)
    a_particle%y = vec(2)
    a_particle%z = vec(3)
  end subroutine

  subroutine move1(a_particle)
    type(particledat), intent(inout) :: a_particle
    type(particledat) :: temp
    call move(a_particle, temp)
    a_particle = temp
  end subroutine

  subroutine move2(oldp, newp)
    type(particledat), intent(in) :: oldp
    type(particledat), intent(out) :: newp
    newp = oldp
    call transmove(oldp%x,oldp%y,oldp%z,newp%x,newp%y,newp%z)
    if (oldp%rod) then
      call rotate(oldp%ux,oldp%uy,oldp%uz,newp%ux,newp%uy,newp%uz)
    end if
  end subroutine
  
  subroutine transmove(xo, yo, zo, xn, yn, zn)
    implicit none
    real(dp), intent(in) :: xo,yo,zo
    real(dp), intent(out) :: xn,yn,zn
    real(dp) :: max_1d 
    max_1d = maxdr/sqrt(3._dp)
    xn = xo + (2.0_dp*grnd()-1.0_dp)*max_1d
    yn = yo + (2.0_dp*grnd()-1.0_dp)*max_1d
    zn = zo + (2.0_dp*grnd()-1.0_dp)*max_1d
  end subroutine

  !Palauttaa partikkelin orientaatiovektorin komponentit
  !sylinterikoordinaatistossa. 
  subroutine unitvec(particle, uro, utheta, uz)
    intrinsic atan2
    type(particledat), intent(in) :: particle
    real(dp), intent(out) :: uro,utheta,uz
    real(dp) :: nx, ny, nz, uxn, uyn, uzn, theta
    theta = -atan2(particle%y, particle%x)
    nx = 0.0_dp
    ny = 0.0_dp
    nz = 1.0_dp
    call xvec2(particle%ux, particle%uy, particle%uz, nx, ny, nz, theta, &
      uxn, uyn, uzn)
    uro = uxn
    utheta = uyn
    uz = uzn
  end subroutine
    
  subroutine rotate(uxo, uyo, uzo, uxn, uyn, uzn)
    real(dp), intent(in) :: uxo,uyo,uzo
    real(dp), intent(out) :: uxn,uyn,uzn
    real(dp) :: theta, nx, ny, nz
    call nvec(nx, ny, nz)
    theta = (2._dp * grnd() - 1._dp) * dthetamax
    call XVEC2(uxo, uyo, uzo, nx, ny, nz, theta, uxn, uyn, uzn)
  end subroutine

  subroutine setmaxmoves(distance, angle)
    implicit none
    real(dp) :: distance, angle
    maxdr = distance/sqrt(3._dp)
    dthetamax = angle
  end subroutine

  subroutine getmaxmoves(distance, angle)
    implicit none
    real(dp), intent(out) :: distance, angle
    distance = maxdr
    angle = dthetamax
  end subroutine 

  pure function max_trans_f(particle) result(mtr)
    type(particledat), intent(in) :: particle
    real(dp) :: mtr
    mtr = maxdr
  end function
  
end module particle


