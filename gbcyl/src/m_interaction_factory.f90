module m_interaction_factory
  use iso_fortran_env
  use class_pair_potential, only: conditional_pair_interaction
  use particle, only: pair_interaction, single_interaction
  use m_gb_interaction, only: gb_interaction
  use m_lj_interaction, only: lj_interaction
  use m_gblj_interaction, only: gblj_interaction
  use json_module, only: json_value, json_get, CK
  use particlewall, only: ljwall_interaction
  implicit none

contains

  function create_pair_interaction(json_val) result(res)
    type(json_value), pointer, intent(in) :: json_val
    class(pair_interaction), pointer :: res
    character(kind=CK, len=:), allocatable :: typestr
    call json_get(json_val, "type", typestr)
    if (typestr == 'conditional_pair_interaction') then
       allocate(res, source=conditional_pair_interaction(json_val))
    else if (typestr == 'gayberne') then
       allocate(res, source=gb_interaction(json_val))
    else if (typestr == 'lj') then
       allocate(res, source=lj_interaction(json_val))
    else if (typestr == 'gblj') then
       allocate(res, source=gblj_interaction(json_val))
    else
       write(error_unit, *) 'ERROR: Unknown interaction type ', typestr
       stop 'create_pair_interaction unable to continue.'
    end if
  end function create_pair_interaction

  function create_single_interaction(json_val) result(res)
    type(json_value), pointer, intent(in) :: json_val
    class(single_interaction), pointer :: res 
    character(kind=CK, len=:), allocatable :: typestr
    res => null()
    call json_get(json_val, "type", typestr)
    if (typestr == 'ljwall_interaction') then
       allocate(res, source=ljwall_interaction(json_val))
    else
       write(error_unit, *) 'ERROR: Unknown interaction type ', typestr
       stop 'create_single_interaction unable to continue.'
    end if
  end function create_single_interaction
  
end module m_interaction_factory
