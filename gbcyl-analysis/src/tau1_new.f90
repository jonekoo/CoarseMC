module tau1
use nrtype, only: dp
use nr, only: mnbrak, brent
implicit none

PRIVATE 


contains

  subroutine tau1(particles, n_particles, value, layer_distance, direction)
    implicit none
    !! 1. choose direction by calculating the order parameter and the 
    !!    corresponding eigenvector
    !! 2. find maximum value of tau1(d) 
  end subroutine tau1



  function tau1(d)
    implicit none
  end function

end module tau1
