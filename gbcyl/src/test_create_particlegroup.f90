program test_create_particlegroup
  use m_particlegroup, only: particlegroup_ptr, particlegroup, &
       create_particlegroup
  use m_rod, only: rod
  use m_particle, only: particle
  use class_poly_box, only: poly_box
  use num_kind, only: dp
  use json_module, only: json_value, json_file, json_initialize
  use m_particlejson_parser, only: particlearray_from_json
  implicit none
  class(particle), allocatable :: particles(:)
  class(particlegroup_ptr), allocatable :: groups(:)
  type(poly_box) :: simbox
  type(json_value), pointer :: json_val
  type(json_file) :: json
  call json_initialize()
  call json%load_file(filename='particlegroup.json')
  call json%get('$', json_val)
  allocate(groups(1))
  call particlearray_from_json(json_val, particles)
  allocate(groups(1)%ptr, source=create_particlegroup(simbox, particles, &
             min_cell_length=6.5_dp, min_boundary_width=0._dp, name='test_group'))
end program test_create_particlegroup
