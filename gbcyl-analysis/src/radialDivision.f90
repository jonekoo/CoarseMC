module division

contains

subroutine radialDivision(radius, div, parray, np, radia, parrays, nps)
use nrtype
use particle
implicit none
real(dp), intent(in) :: radius
integer,intent(in) :: div
type(particledat), dimension(:), intent(in) :: parray
integer, intent(in) :: np
real(dp), dimension(:), pointer :: radia
type(particledat), dimension(:,:), pointer :: parrays
integer, dimension(:), pointer :: nps
real(dp) :: thick,ysqr,xsqr
integer :: ishell, astat,i
integer, dimension(:), allocatable :: endsofarrays
  allocate(radia(div), parrays(div,np), nps(div), endsofarrays(div),&
         & stat=astat)
  if(astat/=0) stop 'radialDivision: memory allocation failed'
!  do i=1,div
!    endsofarrays(i)=1
!  end do
!  thick=radius/div
!  do i=1,np
!    xsqr=parray(i)%x**2
!    ysqr=parray(i)%y**2
!    ishell=ceiling(sqrt(xsqr+ysqr)/thick)
!    parrays(ishell,endsofarrays(ishell))=parray(i)
!  end do
end subroutine radialDivision

end module division
