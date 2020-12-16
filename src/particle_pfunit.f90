!> Unit tests for particle types.
module particle_pfunit
  use m_point, only: point, particlearray_to_json
  use json_module
  use pfunit
  use m_rod, only: rod
  use m_particlejson_parser
  implicit none

contains
  
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
    class(point), allocatable :: in(:)
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
  
  subroutine test_particlearray_json_io
    type(point) :: out(2)
    class(point), allocatable :: in(:)
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
    
  end subroutine test_particlearray_json_io


  subroutine test_point_downcast
    class(point), allocatable :: somepoint
    class(point), allocatable :: another
    integer :: err
    allocate(somepoint, source=point(x=1.0, y=2.0, z=3.0))
    allocate(another, source=point())
    call another%downcast_assign(somepoint, err)
    call assertEqual(0, err, 'point downcast failed')
    select type (somepoint)
    type is (point)
       select type (another)
       type is (point)
          call assertTrue(somepoint == another, 'assignment did not work')
       class default
          call assertTrue(.false., 'another is not point')
       end select
    class default
       call assertTrue(.false., 'somepoint is not point')
    end select
  end subroutine test_point_downcast
  

  subroutine test_rod_downcast
    class(point), allocatable :: somerod
    class(point), allocatable :: another
    integer :: err
    allocate(somerod, source=rod(x=1.0, ux=1.0, uz=0.0))
    allocate(another, source=rod())
    call another%downcast_assign(somerod, err)
    call assertEqual(0, err, 'rod downcast failed')
    select type (somerod)
    type is (rod)
       select type (another)
       type is (rod)
          call assertTrue(somerod == another, 'assignment did not work')
       class default
          call assertTrue(.false., 'another is not rod')
       end select
    class default
       call assertTrue(.false., 'somerod is not rod')
    end select

  end subroutine test_rod_downcast
  
end module particle_pfunit
