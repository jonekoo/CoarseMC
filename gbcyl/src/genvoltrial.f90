module genvoltrial
!! Could make this module independent of the box if only the scaling would be
!! returned and no calls needed for getting the volume.
use class_poly_box
use nrtype
implicit none

contains

  function genvoltrial_scale(simbox, maxscaling, genstate, scalingtype) result(scaling)
    include 'rng.inc'
    type(poly_box), intent(inout) :: simbox
    real(dp), intent(in) :: maxscaling
    type(rngstate), intent(inout) :: genstate
    character(len=*), intent(in) :: scalingtype
    real(dp), dimension(3) :: scaling
    
    if (len(scalingtype) == 1) then
      scaling = genvoltrial_scale1d(simbox, maxscaling, genstate, scalingtype)
    else if (len(scalingtype) == 2) then
      scaling = genvoltrial_scale2d(simbox, maxscaling, genstate, scalingtype) 
    else if (len(scalingtype) == 3) then
      scaling = genvoltrial_scale3d(simbox, maxscaling, genstate)
    else 
      write(*, *) scalingtype
      stop 'genvoltrial_scale: scalingtype not recognized, stopping!'
    end if
  end function

  function genvoltrial_scale1d(simbox, maxscaling, genstate, axis) result(scaling)
    include 'rng.inc'
    type(poly_box), intent(inout) :: simbox
    real(dp), intent(in) :: maxscaling
    type(rngstate), intent(inout) :: genstate
    character, intent(in) :: axis
    real(dp), dimension(3) :: scaling
    real(dp) :: scaling1d
    real(dp) :: Vn, Vo, dV

    dV = (2._dp * real(rng(genstate), dp) - 1._dp) * maxscaling
    Vo = volume(simbox)
    Vn = Vo + dV
    scaling1d = Vn/Vo
    if (axis == 'z') then
      scaling = (/1._dp, 1._dp, scaling1d/)
      call setz(simbox, getz(simbox)*scaling1d)
    else if (axis == 'y') then
      scaling = (/1._dp, scaling1d, 1._dp/)
      call sety(simbox, gety(simbox)*scaling1d)
    else if (axis == 'x') then
      scaling = (/scaling1d, 1._dp, 1._dp/)
      call setx(simbox, getx(simbox)*scaling1d)
    end if
  end function

  function genvoltrial_scale2d(simbox, maxscaling, genstate, plane) result(scaling)
    include 'rng.inc'
    type(poly_box), intent(inout) :: simbox
    real(dp), intent(in) :: maxscaling
    type(rngstate), intent(inout) :: genstate
    character(len=2), intent(in) :: plane
    real(dp), dimension(3) :: scaling
    real(dp) :: scaling1d
    real(dp) :: Vn, Vo, dV

    !! This assumes a rectangular simbox !!
    dV = (2._dp * real(rng(genstate), dp) - 1._dp) * maxscaling
    Vo = volume(simbox)
    Vn = Vo + dV
    scaling1d = sqrt(Vn/Vo)
    if (plane == 'xy' .or. plane == 'yx') then
      call setx(simbox, getx(simbox)*scaling1d)
      call sety(simbox, gety(simbox)*scaling1d)
      scaling = (/scaling1d, scaling1d, 1._dp/)
    else if (plane == 'yz' .or. plane == 'zy') then
      call setz(simbox, getz(simbox)*scaling1d)
      call sety(simbox, gety(simbox)*scaling1d)
      scaling = (/1._dp, scaling1d, scaling1d/)
    else if (plane == 'xz' .or. plane == 'zx') then
      call setz(simbox, getz(simbox)*scaling1d)
      call setx(simbox, getx(simbox)*scaling1d)
      scaling = (/scaling1d, 1._dp, scaling1d/)
    end if
  end function

  function genvoltrial_scale3d(simbox, maxscaling, genstate) result(scaling)
    include 'rng.inc'
    type(poly_box), intent(inout) :: simbox
    real(dp), intent(in) :: maxscaling
    type(rngstate), intent(inout) :: genstate
    real(dp), dimension(3) :: scaling
    real(dp) :: scaling1d
    real(dp) :: Vn, Vo, dV

    !! This assumes a rectangular simbox !!
    dV = (2._dp * real(rng(genstate), dp) - 1._dp) * maxscaling
    Vo = volume(simbox)
    Vn = Vo + dV
    scaling1d = (Vn/Vo)**(1._dp/3._dp)
    call setx(simbox, getx(simbox)*scaling1d)
    call sety(simbox, gety(simbox)*scaling1d)
    call setz(simbox, getz(simbox)*scaling1d)
    scaling = (/scaling1d, scaling1d, scaling1d/)
  end function

end module
