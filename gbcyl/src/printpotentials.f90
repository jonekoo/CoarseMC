program printpotentials
  use cla
  use mc_engine, only: mce_init_json
  use m_nvt_engine, only: simbox, pair_interactions, single_interactions
  use m_point, only: particlearray_wrapper, pair_interaction
  use m_gb_interaction, only: gb_interaction
  use m_lj_interaction, only: lj_interaction
  use m_gblj_interaction, only: gblj_interaction
  use m_lj1wall_interaction, only: lj1wall_interaction
  use m_lj2wall_interaction, only: lj2wall_interaction
  use mpi
  use json_module
  use num_kind, only: dp
  implicit none
  integer, parameter :: id = 0, ntasks=1
  character(len=80) :: input_filename
  character(len=80) :: output_filename
  character(len=80) :: restart_filename
  character(len=80) :: pef_filename
  integer :: err, i, j
  real(dp), parameter :: step = 0.01
  real(dp), parameter :: r(1001) = [(i * step, i = 0, 1000)]
  type(json_value), pointer :: json_val, child_val
  
  call mpi_init(err)
  
  ! Parse command line arguments
  call cla_init
  call cla_register(key='-p', longkey='--pef-file', &
       description='The file to write the potential energy functions.',&
       kkind=cla_char, default='printpotentials.json')
  call cla_register(key='-i', longkey='--input-file', &
       description='input parameter file', &
       kkind=cla_char, default='input-0.json')

  call cla_validate('printpotentials')
  call cla_get('--input-file', input_filename)
  call cla_get('--pef-file', pef_filename)

  ! Initialize
  call mce_init_json(id, ntasks, input_filename, output_filename, &
       restart_filename)
  
  call json_initialize()
  call json_create_object(json_val, '')
  ! compute pair_interactions with sample particles
  do i = 1, size(pair_interactions, 1)
     do j = 1, i
        select type (ia => pair_interactions(i, j)%ptr)
        type is (gb_interaction) 
           call ia%sample(child_val, r)
           call json_add(json_val, child_val)
        type is (lj_interaction)
           call ia%sample(child_val, r)
           call json_add(json_val, child_val)
        type is (gblj_interaction)
           call ia%sample(child_val, r)
           call json_add(json_val, child_val)
        end select
     end do
  end do
  ! compute single_interactions with sample particles
  do i = 1, size(single_interactions)
     select type (ia => single_interactions(i)%ptr)
     type is (lj1wall_interaction)
        call ia%sample(child_val, r, simbox)
        call json_add(json_val, child_val)
     type is (lj2wall_interaction)
        call ia%sample(child_val, r, simbox)
        call json_add(json_val, child_val)
     end select
  end do

  call json_print(json_val, pef_filename)
  call mpi_finalize(err)
  
end program printpotentials
