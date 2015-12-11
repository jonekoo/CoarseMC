module m_interaction_factory
  use class_pair_potential, only: conditional_pair_interaction
  use particle, only: pair_interaction
  use m_gb_interaction, only: gb_interaction
  use json_module, only: json_value, json_get, CK
  implicit none

contains

  function create_pair_interaction(json_val) result(res)
    type(json_value), pointer, intent(in) :: json_val
    class(pair_interaction), pointer :: res
    character(kind=CK, len=:), allocatable :: typestr
    call json_get(json_val, "type", typestr)
    if (typestr == 'conditional_pair_interaction') then
       allocate(res, source=conditional_pair_interaction(json_val))
    else if (typestr == 'gb_interaction') then
       allocate(res, source=gb_interaction(json_val))
    end if
  end function create_pair_interaction
  
end module m_interaction_factory
