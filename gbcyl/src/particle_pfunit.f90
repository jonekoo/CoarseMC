module particle_pfunit
  use particle
  use json_module
  use pfunit
  use m_rod
  use m_point
  use m_particlejson_parser
  implicit none

contains
  
  subroutine test_particledat_json_io
    type(json_value), pointer :: json_val
    type(particledat) :: out, in
    
    out%x = 1.
    out%y = 2.
    out%z = 3.
    out%ux = -1./sqrt(2.)
    out%ux = 1./sqrt(2.)
    out%uz = 0.
    out%rod = .true.
    !! write coordinates to json
    call out%coordinates_to_json(json_val)
    
    !! read coordinates from json
    call in%from_json(json_val)
    call assertTrue(out == in, 'Input and output particledat not equal.')
  end subroutine test_particledat_json_io

  subroutine test_rod_json_io
    type(json_value), pointer :: json_val
    type(rod) :: out, in
    
    out%x = 1.
    out%y = 2.
    out%z = 3.
    out%ux = -1./sqrt(2.)
    out%ux = 1./sqrt(2.)
    out%uz = 0.
    !! write coordinates to json
    call out%coordinates_to_json(json_val)
    !! read coordinates from json
    call in%from_json(json_val)
    call assertTrue(out == in, 'Input and output rod not equal.')
    
  end subroutine test_rod_json_io
  
  subroutine test_point_json_io
    type(json_value), pointer :: json_val
    type(point) :: out, in
    
    out%x = 1.
    out%y = 2.
    out%z = 3.
    !! write coordinates to json
    call out%coordinates_to_json(json_val)
    
    !! read coordinates from json
    call in%from_json(json_val)
    call assertTrue(out == in, 'Input and output point not equal.')    
  end subroutine test_point_json_io

  
  subroutine test_rodarray_json_io
    type(rod) :: out(2)
    class(particledat), allocatable :: in(:)
    type(json_value), pointer :: json_val
    

    call json_create_object(json_val, 'test_particlearray_json_io')
    call particlearray_to_json(json_val, out)
    call particlearray_from_json(json_val, in)
    
    select type(in)
    type is (rod)
       call assertTrue(all(out(:) == in(1:2)), &
            'Input and output particlearrays not equal.')
    class default
       call assertTrue(.false., 'Input is not type rod.')
    end select
    
  end subroutine test_rodarray_json_io
  
  subroutine test_pointarray_json_io
    type(point) :: out(2)
    class(particledat), allocatable :: in(:)
    type(json_value), pointer :: json_val
    

    call json_create_object(json_val, 'test_particlearray_json_io')
    call particlearray_to_json(json_val, out)
    call particlearray_from_json(json_val, in)
    
    select type(in)
    type is (point)
       call assertTrue(all(out(:) == in(1:2)), &
            'Input and output particlearrays not equal.')
    class default
       call assertTrue(.false., 'Input is not of type point.')
    end select
    
  end subroutine test_pointarray_json_io
  
end module particle_pfunit
