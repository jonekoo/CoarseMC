MODULE num_kind
  use iso_fortran_env, only: sp => REAL32, dp => REAL64
  integer, parameter :: spc = kind((1.0_sp, 1.0_sp))
  integer, parameter :: dpc = kind((1.0_dp, 1.0_dp))
END MODULE num_kind
