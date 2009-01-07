!! This code has been used by courtesy of its author Juho Lintuvuori.
!!
!! Jouni Karjalainen, 2008-11-06
!!
!!
subroutine torderparam(N,particle,dmin,dmax,Npoints,latvec,tau,dv)
  use global, only : particledat
  implicit none
  
  integer,intent(in) :: N
  type(particledat),dimension(N),intent(in) :: particle
  double precision,intent(in) :: dmin,dmax
  integer ,intent(in) :: Npoints
  double precision, intent(in) :: latvec(3)
  double precision,intent(out) :: tau,dv

  double precision :: d,deltad,twopi,tau2,arg
  integer :: i,j,nrod
  complex, parameter :: ii=(0.0,1.0)
  complex :: tauhelp

  twopi=8.0*atan(1.0);
  deltad=(dmax-dmin)/real(Npoints);
  nrod=count(particle(1:N)%rod);
  tau=0.0;
  tau2=0.0;
  !write(*,*)'latvec',latvec(1),latvec(2),latvec(3)
  !write(*,*)'dmin,dmax,npoints',dmin,dmax,npoints
  do j=0,Npoints
     d=0.0;
     d = dmin + real(j)*deltad;
     tauhelp=0.0;
     do i=1,N
        if(.not.(particle(i)%rod))cycle;
        arg = latvec(1)*particle(i)%x + latvec(2)*particle(i)%y + &
             latvec(3)*particle(i)%z
        
        tauhelp=tauhelp + exp(ii*arg*twopi/d);
     end do
     !write(*,*)'d,tauhelp',d,tauhelp
     tauhelp=tauhelp/real(nrod)
     tau2=sqrt(tauhelp*conjg(tauhelp));
     if(tau2>tau)then
        tau=tau2;
        dv=d;
     end if
  end do
end subroutine torderparam

subroutine bulkboo(N,Nl,particle,Lx,Ly,Lz,lvec,psi6)
  use global, only : particledat
  implicit none

  integer,intent(in) :: N,Nl
  type(particledat),dimension(N),intent(in) :: particle
  double precision,intent(in) :: Lx,Ly,Lz
  double precision,dimension(3),intent(in) :: lvec
  double precision,intent(out) :: psi6

  integer :: i,j,nrod
  double precision :: wij,sumwij
  complex, parameter :: ii=(0.0,1.0)
  complex :: arg,psi6r
  double precision :: theta
  double precision :: dx,dy,dz,rsq,r,pi,maxval
  double precision,dimension(3) :: rij,u
  double precision,dimension(3),parameter :: zvec=[0.0,0.0,1.0]

  pi=4.0*atan(1.0);
  psi6=0.0;
  psi6r=0.0;
  do i=1,Nl
     sumwij=0.0;
     wij=0.0;
     if(.not.(particle(i)%rod))cycle;
     arg=0.0;
     do j=1,Nl
        if(i==j)cycle;
        if(.not.(particle(j)%rod))cycle;

        dx=particle(i)%x-particle(j)%x;
        dy=particle(i)%y-particle(j)%y;
        dz=particle(i)%z-particle(j)%z;

        !Minimukuva periaate
        dx=dx-nint(dx/Lx)*Lx;
        dy=dy-nint(dy/Ly)*Ly;
        dz=dz-nint(dz/Lz)*Lz;
        
        rsq=dx*dx+dy*dy+dz*dz;
        r=sqrt(rsq);

        !Lasketaan projektio
        ! wij=1 jos r<1.4; wij=0 jos r>1.8; 1.4<r<1.8 lineaarinen int.pol.
        if(r>1.8)then
           cycle;
        else if(r<1.4)then
           wij=1.0;
        else
           wij=-(r-1.4)/(1.8-1.4)+1.0;
        end if
        sumwij=sumwij+wij

        !Lasketaan projektio
        rij(1)=dx;
        rij(2)=dy;
        rij(3)=dz;
       
        rij=rij-(dot_product(rij,lvec))*lvec
        !print *,'Rij',rij
        rij=rij/sqrt(dot_product(rij,rij));
        u(2)=1.0;
        u(1)=0.0;u(3)=0.0;
        u=u-(dot_product(u,lvec))*lvec
        !print *,u,dot_product(u,u)
        !vektorin Rij ja x-akselin valinenkulma theta
        theta=acos(dot_product(rij,u));
      !  print *,'theta',theta*180.0/pi
        arg=arg+wij*exp(ii*6.0*theta);
     end do
     if(sumwij>0.0)then
        !print *,wij,sumwij
        psi6r=psi6r+arg/sumwij;
     end if
  end do
  nrod=count(particle(1:Nl)%rod);
  psi6r=psi6r/nrod;
  psi6=real(psi6r)
end subroutine bulkboo

subroutine orientdistfunc(N,particle,lvec,maxbin,hist)
  use global, only : particledat
  implicit none

  integer,intent(in) :: N
  type(particledat),dimension(N),intent(in) :: particle
  double precision,dimension(3),intent(in) :: lvec
  integer,intent(in) :: maxbin
  double precision,dimension(maxbin),intent(out) :: hist
  
  integer :: nrod,i,bin
  double precision,dimension(3) :: u
  double precision :: costheta,delcostheta,maxval

  delcostheta=1.0/real(maxbin)
  hist(1:maxbin)=0.0;
  do i=1,N
     if(.not.(particle(i)%rod))cycle;
     u(1)=particle(i)%ux;
     u(2)=particle(i)%uy;
     u(3)=particle(i)%uz;

     costheta=dot_product(u,lvec);
     costheta=abs(costheta)
     bin=int(costheta/delcostheta)+1
     !if(costheta<0.0)print *,'hep',costheta,bin

     if(bin<=maxbin)then
        hist(bin)=hist(bin)+1;
     end if
  end do

  maxval=0.0;
  !Normalisoidaan 
  do bin=1,maxbin
     if(hist(bin)>maxval)maxval=hist(bin)
  end do
  hist=hist/maxval;
end subroutine orientdistfunc
