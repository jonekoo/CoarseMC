module m_sphere_interaction
  use num_kind
  use class_parameter_writer
  implicit none
  
  type, abstract :: sphere_interaction
   contains
     procedure(sphere_potential), deferred :: potential
     procedure(sphere_force), deferred :: force
     procedure(sphere_writeparameters), deferred :: writeparameters
  end type sphere_interaction

  interface
     pure subroutine sphere_potential(this, r, energy, overlap)
       import sphere_interaction, dp
       class(sphere_interaction), intent(in) :: this
       real(dp), intent(in) :: r
       real(dp), intent(out) :: energy
       logical, optional, intent(out) :: overlap
     end subroutine sphere_potential

     pure function sphere_force(this, rij)
       import sphere_interaction, dp
       class(sphere_interaction), intent(in) :: this
       real(dp), intent(in) :: rij(3)
       real(dp) :: sphere_force(3)
     end function sphere_force

     subroutine sphere_writeparameters(this, writer)
       import sphere_interaction, parameter_writer
       class(sphere_interaction), intent(in) :: this
       type(parameter_writer), intent(inout) :: writer
     end subroutine sphere_writeparameters
  end interface
  

end module m_sphere_interaction
