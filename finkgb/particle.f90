!Partikkelidatan ja siihen kohdistuvien operaatioiden
!m‰‰rittelyt

module particle
  use nrtype
  implicit none
  
  type particledat
     integer :: index
     real(dp) :: x,y,z,ux,uy,uz
     logical :: rod
  end type particledat

  contains

    subroutine pairV(particlei,particlej,potE,overlap)
    implicit none
    type(particledat),intent(in) :: particlei,particlej
    real(dp), intent(out) :: potE
    logical,intent(out) :: overlap
      potE=0  
    end subroutine pairV
  
    function move(oldparticle) result(newparticle)
    implicit none
    type(particledat), intent(in) :: oldparticle
    type(particledat) :: newparticle
    end function move

end module particle





