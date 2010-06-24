module global
  implicit none
  
  type particledat
     double precision :: x,y,z,ux,uy,uz
     logical :: rod
  end type particledat
  
  double precision, parameter :: kappasig=4.4 !sig_e/sig_s
  double precision, parameter :: kappaeps=20.0 !eps_s/eps_e HUOM! Luckhurst et.al J.Chem.Phys, Vol. 110, No. 14
  double precision, parameter :: myy=1.0;
  double precision, parameter :: nyy=1.0
  double precision, parameter :: sigs=1.0;
  double precision, parameter :: sigs2=1.0*1.0/2.0; ! (sig_s/sqrt(2))^2
  double precision, parameter :: sige2=4.4*4.4/2.0; ! (sig_e/sqrt(2))^2

  double precision, parameter :: epsphers=54.6/20.835  !eps0=20.835 meV; Kr LJ-GB parameter side-side
  double precision, parameter :: epsphere=20.36/20.835 !eps0=20.835 meV; Kr LJ-GB parameter end-end
  double precision, parameter :: epspher=14.13/20.835 ! Kr LJ-epsilon
  double precision, parameter :: rsphere=0.602 ! Kr sigma_LJ/sqrt(2)
  double precision, parameter :: rsphere2=rsphere*rsphere;

end module global

program alkutila
  use global
  implicit none

  !Muuttujat alkutilojen luomista varten
  integer :: Nc,M,iref,N,Nsphere,iatom,astat,ios,mode
  integer :: Napu,i,j,k
  double precision :: u1,u2,s,cell,cell2,dens,pi
  type(particledat) :: base(4)
  type(particledat),dimension(:),allocatable :: particle 
  integer,external :: iargc
  double precision :: Lx,Ly,Lz,L,eps,len
  character(len=40) :: buf
  integer,dimension(:),allocatable :: help
 
  ! mode==2
  integer :: npx,npy,npz
  double precision :: a,d,r,rmin,rmax
  double precision :: dx,dy,dz
  
        Npx=132;
        Npy=132;
        Npz=6;
        read(*, *) npx
        read(*, *) npy
        read(*, *) npz
        N=Npx*Npy*Npz        

        a = 1.1; ! lähinaapuri etäisyys
        d = 3.6; ! kerrostenväli
        
        allocate(particle(N),help(N))
        Nsphere=0;
        particle(1:N)%rod=.true.
        iatom=0;
        do k=0,Npz-1
           do j=0,Npy-1
              do i=0,Npx-1
                 iatom=iatom+1;
                 if(mod(k,2)==0)then
                    if(mod(j,2)==0)then
                       particle(iatom)%x=a*i;
                    else
                       particle(iatom)%x=a*(real(i)+0.5);
                    end if
                 else
                    if(mod(j,2)==0)then
                       particle(iatom)%x=a*(real(i)+0.5);
                    else
                       particle(iatom)%x=a*i;
                    end if
                 end if
                 particle(iatom)%y=sqrt(3.0)/2.0*a*j;
                 particle(iatom)%z=d*k;
                 
                 particle(iatom)%ux=0.0;
                 particle(iatom)%uy=0.0;
                 particle(iatom)%uz=1.0;
                 
              end do
           end do
        end do
       
        Lx = maxval(particle(1:N)%x)+a;
        Ly = maxval(particle(1:N)%y)+a;
        Lz = maxval(particle(1:N)%z)+d;
        rmin=6.0;
        rmax=0.0;
        do i=1,N-1
           do j=i+1,N
              r = (particle(i)%x-particle(j)%x)**2 + &
                   (particle(i)%y-particle(j)%y)**2 + &
                   (particle(i)%z-particle(j)%z)**2;
              if(r<rmin)rmin=r
              if(r>rmax)rmax=r
           end do
        end do
        write(*,*)'iatom',iatom,N
        write(*,*)'Lx...',Lx,Ly,Lz,maxval(particle(1:N)%z)
        write(*,*)'maxval',npx*a,maxval(particle(1:N)%x)
        write(*,*)'dens',real(N)/(Lx*Ly*Lz)
        write(*,*)'rmin',sqrt(rmin),sqrt(rmax)
       
!        call writexyz(particle,N,100,'heksa.xyz');
        
     !Keskitetään koordinaatit laatikon keskustaan
     do i=1,N
        particle(i)%x=particle(i)%x-Lx/2.0 + a/2.0;
        particle(i)%y=particle(i)%y-Ly/2.0 + a/2.0;
        particle(i)%z=particle(i)%z-Lz/2.0 + d/2.0;
     end do
     write(*,*)'Lx/2',Lx/2.0,minval(particle(1:N)%x),maxval(particle(1:N)%x)
     write(*,*)'Ly/2',Ly/2.0,minval(particle(1:N)%y),maxval(particle(1:N)%y)
     write(*,*)'Lz/2',Lz/2.0,minval(particle(1:N)%z),maxval(particle(1:N)%z)

     do iatom=1,N
        if(.not. (particle(iatom)%rod))then
           write(*,*)'x,y,z',particle(iatom)%x,particle(iatom)%y,particle(iatom)%z
        end if
     end do
     rmin=10.0
     do i=1,N-1
        if(particle(i)%rod)cycle
        do j=i+1,N
           if(particle(j)%rod)cycle
           dx=(particle(i)%x-particle(j)%x) 
           dy=(particle(i)%y-particle(j)%y)
           dz=(particle(i)%z-particle(j)%z)
           dx=dx-nint(dx/Lx)*Lx
           dy=dy-nint(dy/Ly)*Ly
           dz=dz-nint(dz/Lz)*Lz
           r=sqrt(dx*dx+dy*dy+dz*dz)
           if(r<rmin)rmin=r
        end do
     end do
     write(*,*)'rmin',rmin
     
     where(particle%rod)
        help=1;
     elsewhere
        help=0
     end where
  open(20,file='atoms.in',status='new',form='formatted',iostat=ios)
  if(ios/=0)stop'Virhe tiedoston atoms.in luomisessa!'
  
  write(20,*)N,Lx,Ly,Lz
  write(20,*)particle(1:N)%x,particle(1:N)%y,particle(1:N)%z,particle(1:N)%ux,particle(1:N)%uy,particle(1:N)%uz,help(1:N)
  close(20);
end program alkutila

subroutine writexyz(particle,N,Nsphere,xyzfile)
  use global
  implicit none

  type(particledat), intent(in) :: particle(*)
  integer, intent(in) :: N,Nsphere
  character(len=40), intent(in) :: xyzfile

  integer :: i,astat,j,Napu
  real, parameter :: sigma0=4.5 ! 4.5 Å; M.A. Bates J.of Chem. Phys. Vol. 120, Nr. 1
  real, parameter :: len2 = 1.5 ! Pituuden puolikas
  double precision :: x1,y1,z1,x2,y2,z2,dist

  open(10,file=xyzfile,iostat=astat)
  if(astat/=0)then
     write(*,*)'Virhe tiedoston, ',trim(xyzfile),' avaamisessa!'
     return;
  end if

  dist=len2/sqrt(3.0);

  write(10,'(1I4)')2*N -Nsphere
  write(10,*)
  do j=1,N;
     if(particle(j)%rod)then
        x1=sigma0*particle(j)%x + dist*particle(j)%ux;
        y1=sigma0*particle(j)%y + dist*particle(j)%uy;
        z1=sigma0*particle(j)%z + dist*particle(j)%uz;
     
        x2=sigma0*particle(j)%x - dist*particle(j)%ux;
        y2=sigma0*particle(j)%y - dist*particle(j)%uy;
        z2=sigma0*particle(j)%z - dist*particle(j)%uz;
        
        write(10,*)'C',x1,y1,z1
        write(10,*)'C',x2,y2,z2
     else
        write(10,*)'In',sigma0*particle(j)%x,sigma0*particle(j)%y,&
             sigma0*particle(j)%z
     end if
  end do
  close(10)
end subroutine writexyz
    
