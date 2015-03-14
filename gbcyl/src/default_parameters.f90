!> This is the main program of the ptgbcyl simulation software for
!! anisotropic particles. 
!!
!! @todo create two different simulation programs: with MPI and without
!!       using e.g. conditional compilation.
!!
program default_parameters
  use mc_engine
  use class_parameter_writer
  use iso_fortran_env
  implicit none
  type(parameter_writer) :: writer
  !! Create writer
  writer = new_parameter_writer(pwunit=output_unit)
  call mce_writeparameters(writer)
end program default_parameters
