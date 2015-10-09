module m_particlegroup
  use particle, only: particledat
  use class_simplelist, only: simplelist
  implicit none
  
  type particlegroup
     type(particledat), allocatable :: particles(:)
     type(simplelist) :: sl
  end type particlegroup

end module m_particlegroup
