module global
  use nrtype
  implicit none
  
  real(dp), parameter :: kappasig=4.4_dp !sig_e/sig_s
  real(dp), parameter :: kappaeps=20._dp !eps_s/eps_e HUOM! Luckhurst et.al J.Chem.Phys, Vol. 110, No. 14
  real(dp), parameter :: myy=1._dp
  real(dp), parameter :: nyy=1._dp
  real(dp), parameter :: sigs=1._dp
  real(dp), parameter :: sigs2=1._dp*1._dp/2._dp ! (sig_s/sqrt(2))^2
  real(dp), parameter :: sige2=4.4_dp*4.4_dp/2._dp ! (sig_e/sqrt(2))^2

  real(dp), parameter :: epsphers=54.6_dp/20.835_dp  !eps0=20.835 meV Kr LJ-GB parameter side-side
  real(dp), parameter :: epsphere=20.36_dp/20.835_dp !eps0=20.835 meV Kr LJ-GB parameter end-end
  real(dp), parameter :: epspher=14.13_dp/20.835_dp ! Kr LJ-epsilon
  real(dp), parameter :: rsphere=0.602_dp ! Kr sigma_LJ/sqrt(2)
  real(dp), parameter :: rsphere2=rsphere*rsphere

end module global

program createbox
  use global
  use particle
  use nrtype
  use class_factory
  use class_poly_box
  implicit none

  !Muuttujat alkutilojen luomista varten
  type(particledat),dimension(:),allocatable :: particles 
  integer,external :: iargc
  real(dp) :: Lx,Ly,Lz
  integer :: ios
  integer :: n
  ! mode==2
  integer :: npx,npy,npz
  type(poly_box) :: simbox
  type(factory) :: thefactory
  character(len=200) :: filename
  Npx=13
  Npy=13
  Npz=6
  write(*, *) 'Give number of particles in directions x, y, z separated by a newline.'
  read(*, *) npx
  read(*, *) npy
  read(*, *) npz
  n=npx*npy*npz
  allocate(particles(n))
  call createcrystal(npx, npy, npz, particles, lx, ly, lz)
  write(filename, '(3(A2,I0),A3)') 'nx', npx, 'ny', npy, 'nz', npz, '.gb' 
  open(20,file=trim(adjustl(filename)),status='new',form='formatted',iostat=ios)
  if(ios/=0) then
    write(6, *) 'Virhe tiedoston ', trim(adjustl(filename)), ' luomisessa!'
    stop 
  end if
  simbox = new_box(lx, ly, lz)
  call writestate(thefactory, 20, simbox, particles) 
  close(20)

contains

subroutine createcrystal(npx, npy, npz, particles, lx, ly, lz)
  integer, intent(in) :: npx, npy, npz
  type(particledat), dimension(npx*npy*npz), intent(out) :: particles 
  real(dp), intent(out) :: lx, ly, lz
  integer :: N,Nsphere,iatom
  integer :: i,j,k
  integer,dimension(:),allocatable :: help
  real(dp) :: a,d,r,rmin,rmax
  real(dp) :: dx,dy,dz
  N=Npx*Npy*Npz        

  a = 1.1_dp 
  d = 3.6_dp ! kerrostenväli
        
  allocate(help(N))
  Nsphere=0
  particles(1:N)%rod=.true.
  iatom=0
  do k=0,Npz-1
    do j=0,Npy-1
      do i=0,Npx-1
        iatom=iatom+1
        if(mod(k,2)==0)then
          if(mod(j,2)==0)then
            particles(iatom)%x=a*real(i,dp)
          else
            particles(iatom)%x=a*(real(i,dp)+0.5_dp)
          end if
        else
          if(mod(j,2)==0)then
            particles(iatom)%x=a*(real(i,dp)+0.5_dp)
          else
            particles(iatom)%x=a*real(i,dp)
          end if
        end if
        particles(iatom)%y=sqrt(3._dp)/2._dp*a*real(j, dp)
        particles(iatom)%z=d*real(k,dp)
                 
        particles(iatom)%ux=0._dp
        particles(iatom)%uy=0._dp
        particles(iatom)%uz=1._dp
                 
      end do
    end do
  end do
       
  Lx = maxval(particles(1:N)%x)+a
  Ly = maxval(particles(1:N)%y)+a
  Lz = maxval(particles(1:N)%z)+d
  rmin=6._dp
  rmax=0._dp
  do i=1,N-1
    do j=i+1,N
      r = (particles(i)%x-particles(j)%x)**2 + &
      (particles(i)%y-particles(j)%y)**2 + &
      (particles(i)%z-particles(j)%z)**2
      if(r<rmin)rmin=r
      if(r>rmax)rmax=r
    end do
  end do
  write(*,*)'iatom',iatom,N
  write(*,*)'Lx...',Lx,Ly,Lz,maxval(particles(1:N)%z)
  write(*,*)'maxval',real(npx, dp)*a,maxval(particles(1:N)%x)
  write(*,*)'dens',real(N,dp)/(Lx*Ly*Lz)
  write(*,*)'rmin',sqrt(rmin),sqrt(rmax)
       
        
  write(*,*) 'Centering coordinates.'
  do i=1,N
    particles(i)%x=particles(i)%x-Lx/2._dp + a/2._dp
    particles(i)%y=particles(i)%y-Ly/2._dp + a/2._dp
    particles(i)%z=particles(i)%z-Lz/2._dp + d/2._dp
  end do
  write(*,*)'Lx/2',Lx/2._dp,minval(particles(1:N)%x),maxval(particles(1:N)%x)
  write(*,*)'Ly/2',Ly/2._dp,minval(particles(1:N)%y),maxval(particles(1:N)%y)
  write(*,*)'Lz/2',Lz/2._dp,minval(particles(1:N)%z),maxval(particles(1:N)%z)

  do iatom=1,N
    if(.not. (particles(iatom)%rod))then
      write(*,*)'x,y,z',particles(iatom)%x,particles(iatom)%y,particles(iatom)%z
    end if
  end do
  rmin=10._dp
  do i=1,N-1
    if(particles(i)%rod)cycle
      do j=i+1,N
        if(particles(j)%rod)cycle
        dx=(particles(i)%x-particles(j)%x) 
        dy=(particles(i)%y-particles(j)%y)
        dz=(particles(i)%z-particles(j)%z)
        dx=dx-anint(dx/Lx)*Lx
        dy=dy-anint(dy/Ly)*Ly
        dz=dz-anint(dz/Lz)*Lz
        r=sqrt(dx*dx+dy*dy+dz*dz)
        if(r<rmin)rmin=r
      end do
  end do
  write(*,*)'rmin',rmin
     
  where(particles%rod)
    help=1
  elsewhere
    help=0
  end where
end subroutine

end program    
