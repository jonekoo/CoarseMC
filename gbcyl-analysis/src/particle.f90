module particle
  use nrtype, only: dp
  use utils, only: rotate_vector

  implicit none


  
  type particledat
     real(dp) :: x, y, z, ux, uy, uz
     logical :: rod
  end type particledat



  contains
  

    
  !! Gives the components of the orientation vector in cylindrical coordinates
  !!
  !! :NOTE: It is not possible to call this routine with 
  !! call unitvec(particle, particle%ux, particle%uy, particle%uz)
  !! 
  subroutine unitvec(particle, uro, utheta, uz)
    implicit none
    intrinsic atan2
    type(particledat), intent(in) :: particle
    real(dp), intent(out) :: uro, utheta, uz
    real(dp) :: nx, ny, nz, theta
    theta = atan2(particle%y, particle%x)
    nx = 0.0
    ny = 0.0 
    nz = 1.0
    call rotate_vector(particle%ux, particle%uy, particle%uz, nx, ny, nz, & 
      & theta, uro, utheta, uz)
  end subroutine unitvec
  


end module particle





