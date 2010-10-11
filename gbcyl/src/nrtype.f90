!! Muokannut kirjasta Numerical Recipes in Fortran 90/95
!! Jouni Karjalainen
!!
module nrtype
  implicit none
  integer, parameter :: I4B=selected_int_kind(9)


  integer, parameter :: sp=selected_real_kind(12)
  integer, parameter :: dp=selected_real_kind(12)

  integer, parameter :: spc=selected_real_kind(12)
  integer, parameter :: dpc=selected_real_kind(12)
  
  integer, parameter :: lgt=kind(.true.)

  

end module nrtype




