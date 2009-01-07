function tau1(rs, n, d) result(t1)
use nrtype
real(dp), dimension(:), intent(in) :: rs
integer, intent(in) :: n
real(dp), intent(in) :: d
real(dp) :: pi
real(dp) :: t1
complex(dpc) :: t1c
complex(dpc), parameter :: i=(0.0, 1.0)
  pi=4*atan(1.0_dp)
  t1c=(0.0, 0.0)
  do k=1,n
    t1c=t1c+exp(2*pi*i*rs(k)/d)
  end do
  t1c=t1c/n
  t1=abs(t1c)
end function tau1
