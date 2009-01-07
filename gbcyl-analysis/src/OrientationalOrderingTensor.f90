module OrientationalOrderingTensor
  
  PRIVATE

  interface eigens
    module procedure eigens_dv, eigens_s
  end interface eigens

  PUBLIC :: eigens

  

  contains 



  subroutine eigens_dv(particlearray,np,d,v)
    use nrtype
    use particle
    implicit none

    type(particledat), dimension(:) :: particlearray
    integer :: np,nrot
    real(dp),dimension(3,3) :: S
    real(dp),dimension(3,3),intent(out) :: v
    real(dp),dimension(3),intent(out) :: d

    call orienmatrix(particlearray,np,S)
    call jacobi(S,3,3,d,v,nrot)    
  end subroutine eigens_dv



  subroutine eigens_s(particleArray, np, eigensString)
    use nrtype
    use particle 
    implicit none

    type(particledat), dimension(:), pointer :: particleArray
    integer, intent(in) :: np
    character(len = *), intent(out) :: eigensString
    real(dp), dimension(3) :: values
    real(dp), dimension(3,3) :: vectors

    call eigens_dv(particleArray, np, values, vectors)    
    write(eigensString, *) values, vectors
  end subroutine eigens_s



  subroutine cylindricalMatrix(particles, N, S)
    use nrtype
    use particle, only : particledat, unitvec
    implicit none

    type(particledat), intent(in) :: particles(*)
    real(dp), intent(out) :: S(3,3)
    real(dp) :: ur,uphi,uz
    integer, intent(in) ::  N
    real(dp) :: Sxx,Sxy,Sxz
    real(dp) :: Syx,Syy,Syz
    real(dp) :: Szx,Szy,Szz
    integer :: i,Nrod
  
    Sxx=0;Sxy=0;Sxz=0;
    Syx=0;Syy=0;Syz=0;
    Szx=0;Szy=0;Szz=0;
    Nrod=0;
    do i=1,N
      if(particles(i)%rod)then
        call unitvec(particles(i),ur,uphi,uz)
        Sxx = Sxx + 3*ur*ur - 1;
        Sxy = Sxy + 3*ur*uphi;
        Sxz = Sxz + 3*ur*uz;
  
        Syx = Syx + 3*uphi*ur;
        Syy = Syy + 3*uphi*uphi - 1;
        Syz = Syz + 3*uphi*uz;
      
        Szx = Szx + 3*uz*ur;
        Szy = Szy + 3*uz*uphi;
        Szz = Szz + 3*uz*uz - 1;

        Nrod=Nrod+1;
      end if
    end do
  
    Sxx=Sxx/(2.0*Nrod); Sxy=Sxy/(2.0*Nrod); Sxz=Sxz/(2.0*Nrod);
    Syx=Syx/(2.0*Nrod); Syy=Syy/(2.0*Nrod); Syz=Syz/(2.0*Nrod);
    Szx=Szx/(2.0*Nrod); Szy=Szy/(2.0*Nrod); Szz=Szz/(2.0*Nrod);

    S(1,1)=Sxx; S(1,2)=Sxy; S(1,3)=Sxz;
    S(2,1)=Syx; S(2,2)=Syy; S(2,3)=Syz;
    S(3,1)=Szx; S(3,2)=Szy; S(3,3)=Szz;
  end subroutine cylindricalMatrix



  subroutine jacobi(a,n,np,d,v,nrot)
  implicit none
  !Numerical recipes 2nd Edition, Volume 1, page 460
  integer :: n,np,nrot ! np>=n
  double precision :: a(np,np),d(np),v(np,np)
  integer, parameter :: nmax=500
  integer :: i,ip,iq,j
  double precision :: c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)

  !Alustetaan v yksikkömatriisiksi
  do ip=1,n
     do iq=1,n
        v(ip,iq)=0.0;
     end do
     v(ip,ip)=1.0;
  end do
  !Alustetaan b ja d a:n diagonaaliksi
  do ip=1,n
     b(ip)=a(ip,ip);
     d(ip)=b(ip);
     z(ip)=0.0;
  end do
  nrot=0;
  i=0;
  do
     i=i+1;
     if(i>=50)then
        write(*,*)'liian monta iteraatiota jacobissa!'
        return;
     end if
     sm=0.0;
     do ip=1,n-1
        do iq=ip+1,n
           sm=sm+abs(a(ip,iq)); !ei diagonaalisten elementtien summa
        end do
     end do
     if(sm==0.0)return;
     if(i < 4)then
        tresh=0.2*sm/(n*n); !Ensimmäiset kolme kierrosta
     else
        tresh=0.0;
     end if
     do ip=1,n-1
        do iq=ip+1,n
           g=100*abs(a(ip,iq));
           !! jos ei diagonaalinen elementti on riittävän pieni 
           !! jätetään kierto väliin neljännen kierroksen jälkeen
           if((i>4) .and. (abs(d(ip))+g==abs(d(ip))) .and. &
               (abs(d(iq))+g==abs(d(iq))))then
              a(ip,iq)=0.0;
           else if(abs(a(ip,iq))>tresh)then
              h=d(iq)-d(ip);  
              if(abs(h)+g == abs(h))then
                 t=a(ip,iq)/h;
              else
                 theta=0.5*h/a(ip,iq);
                 t=1/(abs(theta)+sqrt(1.0+theta*theta));
                 if(theta<0)t=-t;
              end if
              c=1/sqrt(1.0+t*t);
              s=t*c;
              tau=s/(1.0+c);
              h=t*a(ip,iq);
              z(ip)=z(ip)-h;
              z(iq)=z(iq)+h;
              d(ip)=d(ip)-h;
              d(iq)=d(iq)+h;
              a(ip,iq)=0.0;
              do j=1,ip-1
                 g=a(j,ip);
                 h=a(j,iq);
                 a(j,ip)=g-s*(h+g*tau);
                 a(j,iq)=h+s*(g-h*tau);
              end do
              do j=ip+1,iq-1
                 g=a(ip,j);
                 h=a(j,iq);
                 a(ip,j)=g-s*(h+g*tau);
                 a(j,iq)=h+s*(g-h*tau);
              end do
              do j=iq+1,n
                 g=a(ip,j);
                 h=a(iq,j);
                 a(ip,j)=g-s*(h+g*tau);
                 a(iq,j)=h+s*(g-h*tau);
              end do
              do j=1,n
                 g=v(j,ip);
                 h=v(j,iq);
                 v(j,ip)=g-s*(h+g*tau);
                 v(j,iq)=h+s*(g-h*tau);
              end do
              nrot=nrot+1;
           end if
        end do
     end do
     do ip=1,n
        b(ip)=b(ip)+z(ip);
        d(ip)=b(ip);
        z(ip)=0.0;
     end do
  end do
end subroutine jacobi


end module OrientationalOrderingTensor
