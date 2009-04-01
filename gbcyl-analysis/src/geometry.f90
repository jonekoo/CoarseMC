module geometry
use cylinder, only: minimum_image_cyl=>minimum_image

!  interface 
!    function minimum_image(ri, rj, Lx, Ly, Lz) result(rij)
!      use nrtype, only: dp
!      implicit none
!      real(dp), dimension(3) :: rij
!      real(dp), dimension(3), intent(in) :: ri
!      real(dp), dimension(3), intent(in) :: rj
!      real(dp), intent(in) :: Lx, Ly, Lz
!    end function
!  end interface


contains
  
  function minimum_image(ri, rj, Lx, Ly, Lz) result(rij)
    use nrtype, only: dp
    implicit none
    real(dp), dimension(3) :: rij
    real(dp), dimension(3), intent(in) :: ri
    real(dp), dimension(3), intent(in) :: rj
    real(dp), intent(in) :: Lx, Ly, Lz
    rij = minimum_image_cyl(ri, rj, Lx, Ly, Lz)
  end function
    

end module geometry
