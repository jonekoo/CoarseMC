module m_rodsphere_potential
  use num_kind
  use class_parameter_writer
  use json_module
  implicit none
  
  type, abstract :: rodsphere_potential
   contains
     procedure(potential_), deferred :: potential
     procedure(force_), deferred :: force
     procedure(writeparameters_), deferred :: writeparameters
     procedure(to_json_), deferred :: to_json
  end type rodsphere_potential

  abstract interface
     pure subroutine potential_(this, ui, rij, energy, overlap)
       import rodsphere_potential, dp
       class(rodsphere_potential), intent(in) :: this
       real(dp), intent(in) :: ui(3), rij(3)
       real(dp), intent(out) :: energy
       logical, intent(out) :: overlap
     end subroutine potential_

     pure function force_(this, ui, rij)
       import rodsphere_potential, dp
       class(rodsphere_potential), intent(in) :: this
       real(dp), intent(in) :: ui(3), rij(3)
       real(dp) :: force_(3)
     end function force_

     subroutine writeparameters_(this, writer)
       import rodsphere_potential, parameter_writer
       class(rodsphere_potential), intent(in) :: this
       class(parameter_writer), intent(inout) :: writer 
     end subroutine writeparameters_

     subroutine to_json_(this, json_val)
       import rodsphere_potential, json_value
       class(rodsphere_potential), intent(in) :: this
       type(json_value), pointer, intent(inout) :: json_val
     end subroutine to_json_
  end interface
  

end module m_rodsphere_potential
