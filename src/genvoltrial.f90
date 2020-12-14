!> Implements various ways to create a trial volume change in NPT-MC 
!! simulation.
!!
!! @see J. Karjalainen, J. Lintuvuori, V.-V. Telkki, P. Lantto, and 
!!      J. Vaara. Phys. Chem. Chem. Phys., 15:14047, 2013.
!! @see E. de Miguel, E. del Rio, and F. Blas. J. Chem. Phys., 
!!      121(22):11183, 2004.
!!
module genvoltrial
use class_poly_box, only: poly_box, getx, gety, getz, setx, sety, setz
use num_kind, only: dp
use mt_stream, only: rngstate => mt_state, rng => genrand_double1_s
implicit none

contains

!> Scales the simulation box @p simbox in the dimensions specified by
!! @p scalingtype. @p maxscaling sets the maximum volume change.
!! @p genstate holds the random number generator state. Returns the
!! vector describing the scaling of coordinates with (1, 1, 1) being
!! the vector before the scaling.
function genvoltrial_scale(simbox, maxscaling, genstate, scalingtype) &
     result(scaling)
  type(poly_box), intent(inout) :: simbox
  real(dp), intent(in) :: maxscaling
  type(rngstate), intent(inout) :: genstate
  character(len=*), intent(in) :: scalingtype
  real(dp), dimension(3) :: scaling
  
  if (len_trim(adjustl(scalingtype)) == 1) then
     scaling = genvoltrial_scale1d(simbox, maxscaling, genstate, &
          trim(adjustl(scalingtype)))
  else if (len_trim(adjustl(scalingtype)) == 2) then
     scaling = genvoltrial_scale2d(simbox, maxscaling, genstate, &
          trim(adjustl(scalingtype))) 
  else if (len_trim(adjustl(scalingtype)) == 3) then
     scaling = genvoltrial_scale3d(simbox, maxscaling, genstate)
  else 
     write(*, *) trim(adjustl(scalingtype))
     stop 'genvoltrial_scale: scalingtype not recognized, stopping!'
  end if
end function genvoltrial_scale

!> Scales the simulation box @p simbox in a given direction @p axis.
!! @p maxscaling sets the maximum volume change. @p genstate holds the
!! random number generator state.  Returns the vector describing the
!! scaling of coordinates with (1, 1, 1) being the vector before the
!! scaling.
function genvoltrial_scale1d(simbox, maxscaling, genstate, axis) &
     result(scaling)
  type(poly_box), intent(inout) :: simbox
  real(dp), intent(in) :: maxscaling
  type(rngstate), intent(inout) :: genstate
  character, intent(in) :: axis
  real(dp), dimension(3) :: scaling
  real(dp) :: scaling1d
  real(dp) :: Vn, Vo, dV
  real(dp) :: r
  
  call rng(genstate, r)
  dV = (2._dp * real(r, dp) - 1._dp) * maxscaling
  Vo = simbox%volume()
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
end function genvoltrial_scale1d


!> Scales the simulation box @p simbox in the given @p plane. 
!! @p maxscaling sets the maximum volume change. @p genstate holds the
!! random number generator state. Returns the vector describing the
!! scaling of coordinates with (1, 1, 1) being the vector before the
!! scaling.
function genvoltrial_scale2d(simbox, maxscaling, genstate, plane) &
     result(scaling)
  type(poly_box), intent(inout) :: simbox
  real(dp), intent(in) :: maxscaling
  type(rngstate), intent(inout) :: genstate
  character(len=2), intent(in) :: plane
  real(dp), dimension(3) :: scaling
  real(dp) :: scaling1d
  real(dp) :: Vn, Vo, dV
  real(dp) :: r
  
  !! This assumes a rectangular simbox !!
  call rng(genstate, r)
  dV = (2._dp * real(r, dp) - 1._dp) * maxscaling
  Vo = simbox%volume()
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
end function genvoltrial_scale2d


!>  Scales the simulation box @p simbox isotropically in all directions
!! x, y, and z. Returns the vector describing the scaling of
!! coordinates with (1, 1, 1) being the vector before the scaling.
!! @p maxscaling sets the maximum volume change. @p genstate holds the
!! random number generator state. 
function genvoltrial_scale3d(simbox, maxscaling, genstate) result(scaling)
  type(poly_box), intent(inout) :: simbox
  real(dp), intent(in) :: maxscaling
  type(rngstate), intent(inout) :: genstate
  real(dp), dimension(3) :: scaling
  real(dp) :: scaling1d
  real(dp) :: Vn, Vo, dV
  real(dp) :: r
  
  !! This assumes a rectangular simbox !!
  call rng(genstate, r)
  dV = (2._dp * real(r, dp) - 1._dp) * maxscaling
  Vo = simbox%volume()
  Vn = Vo + dV
  scaling1d = (Vn/Vo)**(1._dp/3._dp)
  call setx(simbox, getx(simbox)*scaling1d)
  call sety(simbox, gety(simbox)*scaling1d)
  call setz(simbox, getz(simbox)*scaling1d)
  scaling = (/scaling1d, scaling1d, scaling1d/)
end function genvoltrial_scale3d

end module genvoltrial
