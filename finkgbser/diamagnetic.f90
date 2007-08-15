function diamagnetic(particlei)
  use particle
  use nrtype
  implicit none
  intrinsic dot_product
  type(particledat), intent(in) :: particlei
  real(dp) :: diamagnetic
  real(dp), dimension(3) :: B0,u,Bmol
  real(dp), dimension(3,3) :: sus
  real(dp) :: susiso,anisosus,antipar,par
  real(sp), parameter :: kB=1.38065e-23
  real(sp), parameter :: temp=107.0
  B0=(/0.0_dp,0.0_dp,11.47_dp/)
  susiso=-3290e-30
  anisosus=1750e-30
  antipar=susiso-anisosus/3.0_dp
  par=susiso+2.0_dp*anisosus/3.0_dp
  sus(1:3,1:3)=0
  sus(1,1)=antipar
  sus(2,2)=antipar
  sus(3,3)=par
  u=(/particlei%ux,particlei%uy,particlei%uz/)
  Bmol(3)=dot_product(B0,u)
  Bmol(1)=sqrt(0.5_dp*(dot_product(B0,B0)-Bmol(3)**2))
  Bmol(2)=Bmol(1)
  write (*,*) 'Bmol:',Bmol
  diamagnetic=-dot_product(Bmol,matmul(sus,Bmol))
  diamagnetic=diamagnetic/(kB*temp)
end function
