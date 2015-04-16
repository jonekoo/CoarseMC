module m_rod_interaction
  use nrtype
  use class_parameter_writer
  implicit none
  
  type, abstract :: rod_interaction
   contains
     procedure(rod_potential), deferred :: potential
     procedure(rod_force), deferred :: force
     procedure(rod_writeparameters), deferred :: writeparameters
  end type rod_interaction

  interface
     pure subroutine rod_potential(this, ui, uj, rij, energy, overlap)
       import rod_interaction, dp
       class(rod_interaction), intent(in) :: this
       real(dp), intent(in) :: ui(3), uj(3), rij(3)
       real(dp), intent(out) :: energy
       logical, intent(out) :: overlap
     end subroutine rod_potential
  end interface

  interface
     pure function rod_force(this, ui, uj, rij)
       import rod_interaction, dp
       class(rod_interaction), intent(in) :: this
       real(dp), intent(in) :: ui(3), uj(3), rij(3)
       real(dp) :: rod_force(3)
     end function rod_force
  end interface

  interface
     subroutine rod_writeparameters(this, writer)
       import rod_interaction, parameter_writer
       class(rod_interaction), intent(in) :: this
       type(parameter_writer), intent(inout) :: writer 
     end subroutine rod_writeparameters
  end interface
  

end module m_rod_interaction
