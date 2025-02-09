!> This module provides factory procedures for creating interactions.
module m_interaction_factory
  use iso_fortran_env
  use m_point, only: pair_interaction, single_interaction
  use m_gb_interaction, only: gb_interaction
  use m_lj_interaction, only: lj_interaction
  use m_gblj_interaction, only: gblj_interaction
  use json_module, only: json_value, json_get, CK
  use m_lj1wall_interaction, only: lj1wall_interaction
  use m_lj2wall_interaction, only: lj2wall_interaction
  implicit none

contains

  !> Deserializes a pair interaction from JSON value @p json_val by
  !! redirecting the call to the right constructor.
  !!
  !! @param json_val contains the JSON presentation of the interaction.
  !!
  !! @return pointer to the created pair_interaction.
  !!
  function create_pair_interaction(json_val) result(res)
    type(json_value), pointer, intent(in) :: json_val
    class(pair_interaction), pointer :: res
    character(kind=CK, len=:), allocatable :: typestr
    res => null()
    call json_get(json_val, "type", typestr)
    if (typestr == 'gayberne') then
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

  
  !> Deserializes a single_interaction from json.
  !!
  !! @param json_val contains the JSON presentation.
  !!
  !! @return pointer to the single_interaction created.
  !!
  function create_single_interaction(json_val) result(res)
    type(json_value), pointer, intent(in) :: json_val
    class(single_interaction), pointer :: res 
    character(kind=CK, len=:), allocatable :: typestr
    res => null()
    call json_get(json_val, "type", typestr)
    select case (typestr)
    case ('lj1wall_interaction')
       allocate(res, source=lj1wall_interaction(json_val))
    case ('lj2wall_interaction')
       allocate(res, source=lj2wall_interaction(json_val))
    case default
       write(error_unit, *) 'ERROR: Unknown interaction type ', typestr
       stop 'create_single_interaction unable to continue.'
    end select
  end function create_single_interaction
  
end module m_interaction_factory
