module orientational_ordering
  use nrtype, only: dp
  use particle, only: particledat
  implicit none

  public :: eigens
  public :: orientation_parameter
  public :: resulttype
  public :: to_string
  public :: director

  type resulttype
    real(dp) :: p2
    real(dp), dimension(3) :: director  
  end type resulttype
  
  PRIVATE

  interface eigens
    module procedure eigens_dv
  end interface eigens

  interface orientation_parameter
    module procedure orientation_parameter, orientation_parameter_f
  end interface orientation_parameter
  
  

  contains 


  function director(orientation)
  implicit none
  real(dp), dimension(3) :: director
  type(resulttype) :: orientation
    director = orientation%director
  end function director


  function to_string(res)
  implicit none
  character(len=500) :: to_string 
  type(resulttype), intent(in) :: res
    write(to_string, *) res%p2, res%director
  end function



  function orientation_parameter_f(particles, n_particles)
  implicit none
  type(resulttype) :: orientation_parameter_f
    intrinsic maxloc
    intrinsic maxval
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    integer, dimension(1) :: max_value_position
    integer, parameter :: n_dimensions = 3
    real(dp), dimension(3) :: values
    real(dp), dimension(3, 3) :: vectors
    call eigens_dv(particles, n_particles, values, vectors)
    orientation_parameter_f%p2 = maxval(values)
    max_value_position = maxloc(values)
    orientation_parameter_f%director(1:3) = &
      vectors(1:n_dimensions, max_value_position(1))
  end function orientation_parameter_f



  !! Returns the largest eigenvalue of the orientational ordering tensor and
  !! the eigenvector corresponding to that value. 
  !!
  !! pre-conditions: 
  !! 1. @p particles has @p n_particles > 0 particles.
  !! 
  !! @p particles the particles for which the orientational ordering tensor 
  !! is calculated.
  !! @p n_particles the number of particles.
  !! @p value to be assigned as the orientation parameter.
  !! @p director the eigenvector corresponding to @p value. 
  !!
  subroutine orientation_parameter(particles, n_particles, value, director)
    implicit none
    intrinsic maxloc
    intrinsic maxval
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp), intent(out) :: value
    real(dp), dimension(3), intent(out) :: director
    integer, dimension(1) :: max_value_position
    integer, parameter :: n_dimensions = 3
    real(dp), dimension(3) :: values
    real(dp), dimension(3, 3) :: vectors
    call eigens_dv(particles, n_particles, values, vectors)
    value = maxval(values)
    max_value_position = maxloc(values)
    director(1:n_dimensions) = vectors(1:n_dimensions, max_value_position(1))
  end subroutine orientation_parameter



  subroutine eigens_dv(particles, n_particles, values, vectors)
    implicit none
    type(particledat), dimension(:) :: particles
    integer, intent(in) :: n_particles
    integer :: nrot
    real(dp),dimension(3,3), intent(out) :: vectors
    real(dp),dimension(3), intent(out) :: values
    real(dp),dimension(3,3) :: tensor
    call orientation_tensor(particles, n_particles, tensor)
    call jacobi(tensor, 3, 3, values, vectors, nrot)    
  end subroutine eigens_dv



  subroutine orientation_tensor(particles, n_particles, tensor)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) ::  n_particles
    real(dp), dimension(3, 3), intent(out) :: tensor
    real(dp) :: Sxx, Sxy, Sxz
    real(dp) :: Syx, Syy, Syz
    real(dp) :: Szx, Szy, Szz
    integer :: i, n_rods 
    Sxx = 0; Sxy = 0; Sxz = 0;
    Syx = 0; Syy = 0; Syz = 0;
    Szx = 0; Szy = 0; Szz = 0;
    n_rods = 0;
    do i = 1, n_particles
      if (particles(i)%rod) then
        Sxx = Sxx + 3*particles(i)%ux*particles(i)%ux - 1;
        Sxy = Sxy + 3*particles(i)%ux*particles(i)%uy;
        Sxz = Sxz + 3*particles(i)%ux*particles(i)%uz;
  
        Syx = Syx + 3*particles(i)%uy*particles(i)%ux;
        Syy = Syy + 3*particles(i)%uy*particles(i)%uy - 1;
        Syz = Syz + 3*particles(i)%uy*particles(i)%uz;
      
        Szx = Szx + 3*particles(i)%uz*particles(i)%ux;
        Szy = Szy + 3*particles(i)%uz*particles(i)%uy;
        Szz = Szz + 3*particles(i)%uz*particles(i)%uz - 1;

        n_rods = n_rods+1;
      end if
    end do
  
    Sxx=Sxx/(2.0*n_rods); Sxy=Sxy/(2.0*n_rods); Sxz=Sxz/(2.0*n_rods);
    Syx=Syx/(2.0*n_rods); Syy=Syy/(2.0*n_rods); Syz=Syz/(2.0*n_rods);
    Szx=Szx/(2.0*n_rods); Szy=Szy/(2.0*n_rods); Szz=Szz/(2.0*n_rods);

    tensor(1,1)=Sxx; tensor(1,2)=Sxy; tensor(1,3)=Sxz;
    tensor(2,1)=Syx; tensor(2,2)=Syy; tensor(2,3)=Syz;
    tensor(3,1)=Szx; tensor(3,2)=Szy; tensor(3,3)=Szz;
  end subroutine orientation_tensor



  subroutine jacobi(a, n, np, d, v, nrot)
  implicit none
  !Numerical recipes 2nd Edition, Volume 1, page 460
  integer :: n,np,nrot ! np>=n
  double precision :: a(np,np), d(np), v(np,np)
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
           !! jos ei-diagonaalinen elementti on riittävän pieni 
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



end module orientational_ordering
