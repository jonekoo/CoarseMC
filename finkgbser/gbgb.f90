MODULE gbgb
  use nrtype
  implicit none

  logical, private :: ovrlp=.false.
  real(dp), parameter, private :: kappasig=4.4 !sig_e/sig_s
  real(dp), parameter, private :: kappaeps=20.0 !eps_s/eps_e HUOM! Luckhurst &
  !et.al J.Chem.Phys, Vol. 110, No. 14
  real(dp), parameter, private :: myy=1.0
  real(dp), parameter, private :: nyy=1.0
  real(dp), parameter, private :: sigs=1.0
  real(dp), private, save :: khiie,khiis,khiissq


  contains

  subroutine gbgbV(rij,ui,uj,gbV,ovrlp)
  use nrtype
  implicit none
  real(dp),dimension(3),intent(in) :: rij,ui,uj
  real(dp), intent(out) :: gbV
  logical, pointer :: ovrlp
  real(dp) :: dist,idj,rgb,eps,ids,jds,GB12,GB6
  real(dp), dimension(3) :: su
    gbV=0
    dist=sqrt(dot_product(rij,rij))
    idj=dot_product(ui,uj) 
    su=rij/dist
    ids=dot_product(ui,su)
    jds=dot_product(uj,su)
    eps=ep(idj)*epp(ids,jds,idj);
    rgb=dist-sig(ids,jds,idj)+1.0_dp
    if(rgb<= 0.6)then
      ovrlp=.true.
      return;
    end if
    GB6=(1.0/rgb)**6;
    GB12=GB6*GB6;      
    gbV=GB12-GB6 ;
    gbV=4.0*eps*gbV
  end subroutine gbgbV


    !! Calculates the Gay-Berne potential for two particles 
    !! particlei and particlej. If rgb=rij-sig+sig0<=0.6 
    !! routine returns with ovrlp=.true. 
!    subroutine gbgbV(particlei,particlej,V,ovrlp)
!      use particle, only : particledat
!      use cylinder, only : differences
!      implicit none
!      intrinsic dot_product, sqrt
!      type(particledat),intent(in) :: particlei,particlej
!      real(dp) :: eps,GB6=0.0,GB12=0.0,V
!      logical,pointer :: ovrlp      
!      real(dp) :: idj,ids,jds,rgb,dist
!      real(dp) :: dx,dy,dz      
!      real(dp) :: iux,iuy,iuz
!      real(dp) :: jux,juy,juz
!      real(dp) :: sux,suy,suz
!      call differences(particlei,particlej,dx,dy,dz)
!      V=0.0
!      dist=sqrt(dx*dx+dy*dy+dz*dz)
!      iux=particlei%ux
!      iuy=particlei%uy
!      iuz=particlei%uz
!      jux=particlej%ux
!      juy=particlej%uy
!      juz=particlej%uz
!      idj=dot_product((/iux,iuy,iuz/),(/jux,juy,juz/)) 
!      !!=idotj(particlei,particlej)
!      sux=dx/dist
!      suy=dy/dist
!      suz=dz/dist
!      ids=dot_product((/iux,iuy,iuz/),(/sux,suy,suz/))
!      !!call idotsjdots(particlei,particlej,ids,jds)
!      jds=dot_product((/jux,juy,juz/),(/sux,suy,suz/))
!      eps=ep(idj)*epp(ids,jds,idj);
!      rgb=dist-sig(ids,jds,idj)+1.0_dp
!      if(rgb<= 0.6)then
!        ovrlp=.true.
!        return;
!      end if
!      GB6=(1.0/rgb)**6;
!      GB12=GB6*GB6;      
!      V=GB12-GB6 ;
!      V=4.0*eps*V
!    end subroutine gbgbV


    real(dp) pure function sig(ids,jds,idj)
      implicit none
      intrinsic sqrt
      real(dp),intent(in) :: ids,jds,idj
      real(dp) :: idssq,jdssq,idjsq

      idssq=ids*ids
      jdssq=jds*jds
      idjsq=idj*idj
      sig=1-khiis*(idssq+jdssq-2*khiis*ids*jds*idj)/(1-khiissq * idjsq);
      !write(*,*)'sig=',sig
      sig=1.0/sqrt(sig);
      !write(*,*)'sig=',sig
    end function sig

  
    real(dp) pure function ep(idj)
      implicit none
      intrinsic sqrt
      real(dp),intent(in) :: idj
      ep=1.0/sqrt(1-khiissq*idj**2);
    end function ep


  
    !Gay-Berne potentiaalin parametrit
    real(dp) pure function khiieps()
      khiieps=(kappaeps**(1.0/myy)-1.0)/(kappaeps**(1.0/myy)+1.0);
    end function khiieps     

    real(dp) pure function khiisig()
      khiisig=(kappasig*kappasig-1.0)/(kappasig*kappasig+1.0);
    end function khiisig



    real(dp) pure function epp(ids,jds,idj)
      implicit none
      real(dp),intent(in) :: ids,jds,idj
      epp=1.0_dp-khiie*(ids**2+jds**2-2.0_dp*khiie*ids*jds*idj)/&
          (1.0_dp-khiie**2 * idj**2)
    end function epp


    subroutine initgbgb()
      implicit none
      khiie=khiieps()
      khiis=khiisig()
      khiissq=khiis**2
    end subroutine initgbgb
      
end module gbgb
