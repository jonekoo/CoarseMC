!> Routines to read particles from JSON.
module m_particlejson_parser
  use iso_fortran_env, only: error_unit
  use m_particle, only: particle
  use m_point, only: point
  use m_rod, only: rod
  use json_module, only: json_value, json_get_child, json_count, json_get, CK
  use m_json_wrapper, only: get_parameter
  implicit none
  
contains

  !> Deserializes an array of @p particles from JSON value @p json_val.
  !! A factory routine.
  subroutine particlearray_from_json(json_val, particles)
    type(json_value), pointer, intent(in) :: json_val
    class(particle), allocatable, intent(inout) :: particles(:)
    character(kind=CK, len=:), allocatable :: typestr
    integer :: n, i
    type(json_value), pointer :: coordinates_json, particle_json
    call get_parameter(json_val, 'type', typestr)
    !call get_parameter(json_val, 'size', n, error_lb=0)
    call json_get(json_val, 'coordinates', coordinates_json)
    n = json_count(coordinates_json)
    if (typestr == 'rod') then
       allocate(particles(n), source=rod())
    else if (typestr == 'point') then
       allocate(particles(n), source=point())
    else
       write(error_unit, *) 'ERROR: particlearray_from_json: particle type ', &
            typestr, ' not recognized. Stopping.'
    end if
    do i = 1, size(particles)
       call json_get_child(coordinates_json, i, particle_json)
       call particles(i)%from_json(particle_json) 
    end do
  end subroutine particlearray_from_json
  
end module m_particlejson_parser
