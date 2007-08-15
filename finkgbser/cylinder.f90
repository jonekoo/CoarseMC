module cylinder
  use nrtype

  real(dp),private, save :: R
  real(dp), private, save :: Lz,T,P
  
  contains

  !Palauttaa tämän sylinterin tilavuuden
  real(dp) pure function volume()
    implicit none
    intrinsic atan
    real(dp) :: pi
    pi=4.0*atan(1.0)
    volume = pi*R**2*Lz
  end function volume

  !Muuttaa tämän sylinterin sädettä
  subroutine setR(newR)
    implicit none
    real(dp),intent(in) :: newR
    R=newR 
  end subroutine setR


  subroutine setLz(newLz)
    implicit none
    real(dp), intent(in) :: newLz
    Lz=newLz
  end subroutine setLz




  !Skaalaa koordinaatit uudelle R:lle, kun vanha R on R0
  subroutine scale(x0,y0,xnew,ynew,R0,Rnew)
    implicit none
    real(dp),intent(in) :: x0, y0, R0,Rnew
    real(dp),intent(out) :: xnew, ynew
   
    xnew=Rnew/R0*x0
    ynew=Rnew/R0*y0
      
  end subroutine scale



  !Palauttaa etäisyyden sylinterin akselilta
  real(dp) pure function rdist(x,y)
    implicit none
    real(dp), intent(in) :: x,y
    rdist=sqrt(x**2+y**2)
  end function rdist


  !Palauttaa annetut karteesiset koordinaatit
  !sylinterikoordinaateissa. z-koordinaatti on
  !tässä sylinterissä perioidinen
  subroutine cylindrical(x,y,z,r,phi,zc)
    implicit none
    intrinsic atan2,nint
    real(dp), intent(in) :: x,y,z
    real(dp), intent(out) :: r, phi, zc 

    r=rdist(x,y)
    phi=atan2(y,x)
    zc=z-nint(z/Lz)*Lz
  end subroutine cylindrical



  !Sylinterin ominaisuuksiin liittyviä aliohjelmia ja funktioita. 

  !Palauttaa tämän sylinterin säteen
  pure function getradius() result(radius)
    implicit none
    real(dp) :: radius
    radius=R
  end function getradius

  !Palauttaa tämän sylinterin korkeuden
  pure function getHeight() result(height)
    implicit none
    real(dp) :: height
    height=Lz
  end function getHeight

  !Palauttaa lämpötilan sylinterissä
  pure function getTemp() result(Temp)
    implicit none
    real(dp) :: Temp
    Temp=T
  end function getTemp

  !Palauttaa paineen sylinterissä
  pure function getPres() result(Pres)
    implicit none
    real(dp) :: Pres
    Pres=P
  end function getPres

  !Alustaa sylinterin säteen ja korkeuden
  subroutine initcylinder(radius,height,temp,pres)
    implicit none
    real(dp),intent(in) :: radius,height,temp,pres
    R=radius
    Lz=height
    T=temp
    P=pres
  end subroutine initcylinder


  pure function shellVolume(rad,drad) result(vol)
    implicit none
    intrinsic atan
    real(dp),intent(in) :: rad,drad
    real(dp) :: Pi,vol
 
    Pi=4*atan(1.0)
    vol=2.0*Pi*Lz*rad*drad
  end function shellVolume

end module cylinder
