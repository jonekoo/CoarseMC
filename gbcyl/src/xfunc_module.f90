module xfunc_module
  use m_particledat
  use num_kind
  implicit none

contains

  function rho(p) 
  real(dp) :: rho
  type(particledat), intent(in) :: p
    real(dp), dimension(3) :: pos
    pos = p%position()
    rho = sqrt(pos(1)**2 + pos(2)**2)
  end function rho

end module xfunc_module
