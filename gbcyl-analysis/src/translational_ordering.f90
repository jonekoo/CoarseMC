module translational_ordering
  use nrtype, only: dp, dpc
  use particle, only: particledat




  PUBLIC :: tau1

  interface tau1
    module procedure tau1_string, tau1_double
  end interface tau1



  PRIVATE 

  contains



  subroutine tau1_string(particles, n_particles, direction, & 
                       & tau1_s)
    implicit none 
    type(particledat), dimension(:), pointer :: particles
    integer, intent(in) :: n_particles
    real(dp), dimension(3) :: direction
    character(len = 78), intent(out) :: tau1_s
    real(dp) :: tau1_d
    real(dp) :: layer_distance
    call tau1_double(particles, n_particles, direction, tau1_d, layer_distance)
    write(tau1_s, *) tau1_d, layer_distance 
  end subroutine tau1_string



  subroutine tau1_double(particles, n_particles, direction, & 
                         & tau1_d, layer_distance)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp), dimension(3), intent(in) :: direction
    real(dp), intent(out) :: tau1_d
    real(dp), intent(out) :: layer_distance
    real(dp), parameter :: d_layer_min = 2.2    
    real(dp), parameter :: d_layer_max = 7.0
    integer, parameter :: n_points = 100
    call torderparam(n_particles, particles, d_layer_min, d_layer_max, & 
                   & n_points, direction, tau1_d, layer_distance) 
  end subroutine tau1_double




  !! This code has been used by courtesy of its author Juho Lintuvuori.
  !!
  !! Modifications by Jouni Karjalainen, 2008-11-06
  !!
  !!
  subroutine torderparam(N, particle, dmin, dmax, Npoints, latvec, tau, dv)
    implicit none
    integer,intent(in) :: N
    type(particledat),dimension(N),intent(in) :: particle
    real(dp),intent(in) :: dmin, dmax
    integer ,intent(in) :: Npoints
    real(dp), intent(in) :: latvec(3)
    real(dp),intent(out) :: tau, dv
    real(dp) :: d, deltad, twopi, tau2, arg
    integer :: i, j, nrod
    complex, parameter :: ii = (0.0, 1.0)
    complex :: tauhelp
    twopi=8.0*atan(1.0);
    deltad=(dmax-dmin)/real(Npoints);
    nrod=count(particle(1:N)%rod);
    tau=0.0;
    tau2=0.0;
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
      tauhelp=tauhelp/real(nrod)
      tau2=sqrt(tauhelp*conjg(tauhelp));
      if(tau2>tau)then
        tau=tau2;
        dv=d;
      end if
    end do
  end subroutine torderparam



  function tau1_simple(rs, n, d) result(t1)
    implicit none
    real(dp), dimension(:), intent(in) :: rs
    integer, intent(in) :: n
    real(dp), intent(in) :: d
    real(dp) :: pi
    real(dp) :: t1
    complex(dpc) :: t1c
    complex(dpc), parameter :: i=(0.0, 1.0)
    integer :: k
    pi = 4*atan(1.0_dp)
    t1c = (0.0, 0.0)
    do k = 1, n
      t1c = t1c + exp(2*pi*i*rs(k)/d)
    end do
    t1c = t1c/n
    t1 = abs(t1c)
  end function tau1_simple

end module translational_ordering
