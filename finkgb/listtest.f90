program listtest
use linkedlist
use particle
use nrtype
implicit none
type(listtype), pointer :: list
type(particledat),dimension(:), allocatable :: parray
integer :: i 
allocate(parray(5))
do i=1,5
  parray(i)%x=i
  parray(i)%y=0
  parray(i)%z=0
  parray(i)%ux=1
  parray(i)%uy=0
  parray(i)%uz=0
  parray(i)%rod=.true.
  call prependToList(list,parray(i))
end do
call printList(list)

end program listtest
