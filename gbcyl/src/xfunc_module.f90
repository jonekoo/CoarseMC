module xfunc_module
implicit none
contains
  function rho(p) 
  use particle
  use num_kind
  implicit none
  real(dp) :: rho
  type(particledat), intent(in) :: p
    real(dp), dimension(3) :: pos
    pos = position(p)
    rho = sqrt(pos(1)**2 + pos(2)**2)
  end function rho

end module xfunc_module
