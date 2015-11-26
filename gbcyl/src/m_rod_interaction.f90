module m_rod_interaction
  use num_kind
  use class_parameter_writer
  use json_module, only: json_value
  implicit none
  
  type, abstract :: rod_interaction
   contains
     procedure(rod_potential), deferred :: potential
     procedure(rod_force), deferred :: force
     procedure(rod_writeparameters), deferred :: writeparameters
     procedure(to_json), deferred :: to_json
  end type rod_interaction

  abstract interface
     pure subroutine rod_potential(this, ui, uj, rij, energy, overlap)
       import rod_interaction, dp
       class(rod_interaction), intent(in) :: this
       real(dp), intent(in) :: ui(3), uj(3), rij(3)
       real(dp), intent(out) :: energy
       logical, intent(out) :: overlap
     end subroutine rod_potential

     pure function rod_force(this, ui, uj, rij)
       import rod_interaction, dp
       class(rod_interaction), intent(in) :: this
       real(dp), intent(in) :: ui(3), uj(3), rij(3)
       real(dp) :: rod_force(3)
     end function rod_force

     subroutine rod_writeparameters(this, writer)
       import rod_interaction, parameter_writer
       class(rod_interaction), intent(in) :: this
       type(parameter_writer), intent(inout) :: writer 
     end subroutine rod_writeparameters

     subroutine to_json(this, json_val)
       import rod_interaction, json_value
       class(rod_interaction), intent(in) :: this
       type(json_value), pointer, intent(inout) :: json_val
     end subroutine to_json
  end interface
  

end module m_rod_interaction
