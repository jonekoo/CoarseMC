pure function zhangI(m, k, rc)
  use num_kind, only: dp
  use cylinder_integrals, only: hyp2f1
  use m_hyp2f1_negint, only: hyp2f1_negint
  implicit none
  integer, intent(in) :: m
  real(dp), intent(in) :: k, rc
  real(dp) :: zhangI
  real(dp), parameter :: pi = 4._dp * atan(1._dp)
  real(dp) :: h2f1
  if (mod(m, 2) == 0) then
     h2f1 = hyp2f1(-(m - 1) / 2._dp, -(m - 3) / 2._dp, 1._dp, k**2) 
  else
     h2f1 = hyp2f1_negint(-(m - 1) / 2, -(m - 3) / 2, 1, k**2)
  end if
  zhangI = 2 / ((m - 1) * rc**(m - 1) * (1 - k**2)**(m - 1)) * pi * h2f1
end function zhangI

!> Derivative of zhangI with respect to k.
pure function d_zhangI(m, k, rc)
  use num_kind, only: dp
  use cylinder_integrals, only: hyp2f1
  use m_hyp2f1_negint, only: hyp2f1_negint
  implicit none
  integer, intent(in) :: m
  real(dp), intent(in) :: k, rc
  real(dp) :: d_zhangI
  real(dp), parameter :: pi = 4 * atan(1._dp)
  if (mod(m, 2) /= 0 .and. m > 3) then
     d_zhangI = -k * (1 - k**2)**(-m) * pi * rc**(1 - m) * &
          (-4 * hyp2f1_negint((1 - m) / 2, (3 - m) / 2, 1, k**2) + &
          (-1 + k**2) * (-3 + m) * &
          hyp2f1_negint((3 - m) / 2, (5 - m) / 2, 2, k**2))
  else
     d_zhangI = -k * (1 - k**2)**(-m) * pi * rc**(1 - m) * &
          (-4 * hyp2f1((1 - m) / 2._dp, (3 - m) / 2._dp, 1._dp, k**2) + &
          (-1 + k**2) * (-3 + m) * &
          hyp2f1((3 - m) / 2._dp, (5 - m) / 2._dp, 2._dp, k**2))
  end if
end function d_zhangI


pure function repwall2(k, rc)
  use num_kind, only: dp
  use cylinder_integrals, only: zhangI
  implicit none
  real(dp), intent(in) :: k, rc
  real(dp) :: repwall2
  real(dp), parameter :: pi = 4._dp * atan(1._dp)
  repwall2 = 63. * pi / 64 * zhangI(10, k, rc)
end function repwall2


pure function attwall2(k, rc)
  use num_kind, only: dp
  use cylinder_integrals, only: zhangI
  implicit none
  real(dp), intent(in) :: k, rc
  real(dp) :: attwall2
  real(dp), parameter :: pi = 4._dp * atan(1._dp)
  attwall2 = 3. * pi  / 2 * zhangI(4, k, rc)
end function attwall2


pure function d_repwall(k, rc)
  use num_kind, only: dp
  use cylinder_integrals, only: d_zhangI
  implicit none
  real(dp), intent(in) :: k, rc
  real(dp) :: d_repwall
  real(dp), parameter :: pi = 4._dp * atan(1._dp)
  d_repwall = 63. * pi / 64. * d_zhangI(10, k, rc)
end function d_repwall

pure function d_attwall(k, rc)
  use num_kind, only: dp
  use cylinder_integrals, only: d_zhangI
  implicit none
  real(dp), intent(in) :: k, rc
  real(dp) :: d_attwall
  real(dp), parameter :: pi = 4. * atan(1._dp)
  d_attwall = 3. * pi / 2. * d_zhangI(4, k, rc)
end function d_attwall


