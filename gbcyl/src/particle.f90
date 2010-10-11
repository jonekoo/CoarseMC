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
  public :: createparticle
  public :: writeparticle
  public :: position
  public :: orientation
  public :: particle_writeparameters
  public :: getmaxmoves
  public :: move
  public :: setmaxmoves
  public :: setposition
  public :: setorientation
  public :: maxtrans
  
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
  real(dp), save :: maxdr = -1._dp
  character(len = *), parameter :: typeid = 'gb'

  interface initparticle
    module procedure init_particleold, init_particleparameterizer
  end interface

  interface maxtrans
    module procedure maxtransf
  end interface

  interface move
    module procedure move1, move2
  end interface

  contains

  subroutine init_particleparameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'max_translation', maxdr)
    call getparameter(reader, 'max_rotation', dthetamax)
    if(maxdr < 0._dp) stop 'No max_translation given. Stopping.'
    if(dthetamax < 0._dp) stop 'No max_rotation given. Stopping.'
  end subroutine

  subroutine init_particleold(maxTranslation, maxRotation)
    implicit none
    real(dp), intent(in) :: maxTranslation 
    real(dp), intent(in) :: maxRotation
    dthetamax = maxRotation
    maxdr = maxTranslation
  end subroutine

  subroutine particle_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'Particle parameters')
    call writeparameter(writer, 'max_translation', maxdr)
    call writeparameter(writer, 'max_rotation', dthetamax)    
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

  function createparticle(particlestring) result(pp)
    type(particledat) :: pp
    character(len = *), intent(in) :: particlestring
    character(len = 3) :: typeid
    read(particlestring, *) typeid
    if ('gb' == typeid) then
      read(particlestring, *) typeid, pp%x, pp%y, pp%z, pp%ux, pp%uy, pp%uz
      pp%rod = .true.
    else
      stop 'Error: Could not read particle from string. Stopping.'
    end if
  end function

  subroutine writeparticle(aparticle, writeunit)
    integer, intent(in) :: writeunit
    type(particledat), intent(in) :: aparticle
    write(writeunit, '(A, 6' // fmt_char_dp() // ')', advance='no') &
    typeid, aparticle%x, aparticle%y, aparticle%z, aparticle%ux, &
    aparticle%uy, aparticle%uz 
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

  subroutine move1(aparticle)
    type(particledat), intent(inout) :: aparticle
    type(particledat) :: temp
    call move(aparticle, temp)
    aparticle = temp
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
    real(dp), intent(in) :: xo, yo, zo
    real(dp), intent(out) :: xn, yn, zn
    real(dp) :: max1d 
    max1d = maxdr/sqrt(3._dp)
    xn = xo + (2._dp * grnd() - 1._dp) * max1d
    yn = yo + (2._dp * grnd() - 1._dp) * max1d
    zn = zo + (2._dp * grnd() - 1._dp) * max1d
  end subroutine

  !Palauttaa partikkelin orientaatiovektorin komponentit
  !sylinterikoordinaatistossa. 
  !subroutine unitvec(aparticle, uro, utheta, uz)
  !  intrinsic atan2
  !  type(particledat), intent(in) :: aparticle
  !  real(dp), intent(out) :: uro,utheta,uz
  !  real(dp) :: nx, ny, nz, uxn, uyn, uzn, theta
  !  theta = -atan2(aparticle%y, aparticle%x)
  !  nx = 0.0_dp
  !  ny = 0.0_dp
  !  nz = 1.0_dp
  !  call xvec2(aparticle%ux, aparticle%uy, aparticle%uz, nx, ny, nz, theta, &
  !    uxn, uyn, uzn)
  !  uro = uxn
  !  utheta = uyn
  !  uz = uzn
  !end subroutine
    
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
    real(dp), intent(in) :: distance, angle
    maxdr = distance
    dthetamax = angle
  end subroutine

  subroutine getmaxmoves(distance, angle)
    implicit none
    real(dp), intent(out) :: distance, angle
    distance = maxdr
    angle = dthetamax
  end subroutine 

  pure function maxtransf(aparticle) result(mtr)
    type(particledat), intent(in) :: aparticle
    real(dp) :: mtr
    mtr = maxdr
  end function
  
end module
