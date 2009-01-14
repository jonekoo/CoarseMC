module nr

  INTERFACE 
    FUNCTION brent(ax,bx,cx,func,tol,xmin) 
    USE nrtype 
    REAL(SP), INTENT(IN) :: ax,bx,cx,tol 
    REAL(SP), INTENT(OUT) :: xmin 
    REAL(SP) :: brent 
    INTERFACE 
      FUNCTION func(x) 
      USE nrtype 
      REAL(SP), INTENT(IN) :: x 
      REAL(SP) :: func 
      END FUNCTION func 
    END INTERFACE 
    END FUNCTION brent 
  END INTERFACE 

  interface
    subroutine bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
      use nrtype
      real(sp), dimension(:), intent(inout) :: y
      real(sp), dimension(:), intent(in) :: dydx,yscal
      real(sp), intent(inout) :: x
      real(sp), intent(in) :: htry,eps
      real(sp), intent(out) :: hdid,hnext
      interface 
        subroutine derivs(x,y,dydx)
        use nrtype
        real(sp), intent(in) :: x
        real(sp), dimension(:), intent(in) :: y
        real(sp), dimension(:), intent(out) :: dydx
        end subroutine derivs
      end interface
    end subroutine bsstep
  end interface

  interface elle
    function elle_s(phi,ak)
    use nrtype
    real(sp), intent(in) :: phi,ak
    real(sp) :: elle_s
    end function elle_s
  
    function elle_v(phi,ak)
    use nrtype
    real(sp), dimension(:), intent(in) :: phi,ak
    real(sp), dimension(size(phi)) :: elle_v
    end function elle_v
  end interface

  interface ellf
    function ellf_s(phi,ak)
    use nrtype
    real(sp), intent(in) :: phi,ak
    real(sp) :: ellf_s
    end function ellf_s

    function ellf_v(phi,ak)
    use nrtype
    real(sp),dimension(:),intent(in) :: phi,ak
    real(sp),dimension(size(phi)) :: ellf_v
    end function ellf_v
  end interface
  
  interface
    subroutine hypdrv(s,ry,rdyds)
      use nrtype
      real(sp), intent(in) :: s
      real(sp), dimension(:), intent(in) :: ry
      real(sp), dimension(:), intent(out) :: rdyds
    end subroutine hypdrv 
  end interface

  interface 
    function hypgeo(a,b,c,z)
      use nrtype
      complex(spc), intent(in) :: a,b,c,z
      complex(spc) :: hypgeo
    end function hypgeo
  end interface

  interface
    subroutine hypser(a,b,c,z,series,deriv)
      use nrtype
      complex(spc), intent(in) :: a,b,c,z
      complex(spc), intent(out) :: series, deriv
    end subroutine hypser
  end interface

  interface 
    function locate(xx,x)
      USE nrtype
      real(sp), dimension(:), intent(in) :: xx
      real(sp), intent(in) :: x
      integer(I4B) :: locate
    end function locate
  end interface


  INTERFACE
    SUBROUTINE mmid(y,dydx,xs,htot,nstep,yout,derivs)
      USE nrtype
      INTEGER(I4B), INTENT(IN) :: nstep
      REAL(SP), INTENT(IN) :: xs,htot
      REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
      REAL(SP), DIMENSION(:), INTENT(OUT) :: yout

      INTERFACE
        SUBROUTINE derivs(x,y,dydx)
          USE nrtype
          REAL(SP), INTENT(IN) :: x 
          REAL(SP), DIMENSION(:), INTENT(IN) :: y
          REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
        END SUBROUTINE derivs
      END INTERFACE
    END SUBROUTINE mmid
  END INTERFACE

  INTERFACE 
    SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func) 
    USE nrtype 
    REAL(SP), INTENT(INOUT) :: ax,bx 
    REAL(SP), INTENT(OUT) :: cx,fa,fb,fc 
    INTERFACE 
      FUNCTION func(x) 
      USE nrtype 
      REAL(SP), INTENT(IN) :: x 
      REAL(SP) :: func 
      END FUNCTION func 
    END INTERFACE 
    END SUBROUTINE mnbrak 
  END INTERFACE 

  interface 
    subroutine odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
      use nrtype
      real(sp), dimension(:), intent(inout) :: ystart
      real(sp), intent(in) :: x1,x2,eps,h1,hmin
      interface 
        subroutine derivs(x,y,dydx)
          use nrtype
          real(sp), intent(in) :: x
          real(sp), dimension(:), intent(in) :: y
          real(sp), dimension(:), intent(out) :: dydx
        end subroutine derivs
        
        subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
          use nrtype
          real(sp), dimension(:), intent(inout) :: y
          real(sp), dimension(:), intent(in) :: dydx,yscal
          real(sp), intent(inout) :: x
          real(sp), intent(in) :: htry, eps
          real(sp), intent(out) :: hdid, hnext
          interface
            subroutine derivs(x,y,dydx)
              use nrtype
              real(sp), intent(in) :: x
              real(sp), dimension(:), intent(in) :: y 
              real(sp), dimension(:), intent(out) :: dydx
            end subroutine derivs
          end interface
        end subroutine rkqs
      end interface
    end subroutine odeint
  end interface



  INTERFACE
    SUBROUTINE pzextr(iest,xest,yest,yz,dy)
      USE nrtype
      INTEGER(I4B), INTENT(IN) :: iest
      REAL(SP), INTENT(IN) :: xest
      REAL(SP), DIMENSION(:), INTENT(IN) :: yest
      REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
    END SUBROUTINE pzextr
  END INTERFACE

  interface rd
    function rd_s(x,y,z)
    use nrtype
    real(sp), intent(in) :: x,y,z
    real(sp) :: rd_s
    end function rd_s

    function rd_v(x,y,z)
    use nrtype
    real(sp),dimension(:),intent(in) :: x,y,z
    real(sp),dimension(size(x)) :: rd_v
    end function rd_v
  end interface

  interface rf
    function rf_s(x,y,z)
    use nrtype
    real(sp),intent(in) :: x,y,z
    real(sp) :: rf_s
    end function rf_s

    function rf_v(x,y,z)
    use nrtype
    real(sp), dimension(:), intent(in) :: x,y,z
    real(sp), dimension(size(x)) :: rf_v
    end function
  end interface

  interface 
    subroutine rkck(y,dydx,x,h,yout,yerr,derivs)
      use nrtype
      real(sp), dimension(:), intent(in) :: y,dydx
      real(sp), intent(in) :: x,h
      real(sp), dimension(:), intent(out) :: yout,yerr
      
      interface 
        subroutine derivs(x,y,dydx)
          use nrtype
          real(sp), intent(in) :: x
          real(sp), dimension(:), intent(in) :: y
          real(sp), dimension(:), intent(out) :: dydx
        end subroutine derivs
      end interface
    end subroutine rkck
  end interface

  interface
    subroutine spline(x,y,yp1,ypn,y2)
      USE nrtype
      real(sp), dimension(:), intent(in) :: x,y
      real(sp), intent(in) :: yp1,ypn
      real(sp), dimension(:), intent(out) :: y2
    end subroutine spline
  end interface

  interface
    function splint(xa,ya,y2a,x)
    USE nrtype
    real(sp), dimension(:), intent(in) :: xa,ya,y2a
    real(sp), intent(in) :: x
    real(sp) :: splint
    end function splint 
  end interface

  !! On a purely serial machine, for greater efficiency, remove
  !! the generic name tridag from the following interface,
  !! and put it on the next one after that. 
  !! Done here. -JK
  interface tridag 
    subroutine tridag_ser(a,b,c,r,u)
    USE nrtype
    real(sp), dimension(:), intent(in) :: a,b,c,r
    real(sp), dimension(:), intent(out) :: u
    end subroutine tridag_ser
  end interface


end module nr
