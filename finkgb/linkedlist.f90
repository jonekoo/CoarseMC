module linkedlist
use particle
use particlearray, only : allocateParticleArray, reallocateParticleArray

type listtype
  type(particledat) :: prtcl
  type(listtype), pointer :: next
end type listtype

integer, parameter :: cellmax=1000

contains 

!Adds to the beginning of the list
subroutine prependToList(list,prtcl)
implicit none
type(listtype), pointer :: list,newlist,present
type(particledat), intent(in) :: prtcl
integer :: astat
  present=>list
  if(associated(list)) then
    allocate(newlist,stat=astat)
    if(astat/=0) stop 'addToList: allocation failed'
    newlist%prtcl=prtcl
    newlist%next=>present
    list=>newlist
  else
    allocate(list,stat=astat)
    if(astat/=0) stop 'addToList: allocation failed'
    list%prtcl=prtcl
    nullify(list%next)
  end if    
end subroutine prependToList



subroutine freeList(list)
  type(listtype), pointer :: list,next
  do while(associated(list))
    next=>list%next
    deallocate(list)
    nullify(list)
    list=>next
  end do
end subroutine freeList 



subroutine printList(list)
implicit none
type(listtype), pointer :: list,present
  present=>list
  do while(associated(present))
    write(*,*) present%prtcl
    present=>present%next
  end do  
end subroutine printList


subroutine toArray(list,parray,np)
implicit none
type(listtype),pointer :: list,present
type(particledat), dimension(:), pointer :: parray, tmpparray
integer, intent(out) :: np
integer :: k,np0
  np=0
  if(.not. associated(list)) return;
  tmpparray=>allocateParticleArray(500)
  np0=cellmax
  k=0
  present=>list
  do while(associated(present))
    k=k+1
    if(k>np0) then
      tmpparray=>reallocateParticleArray(tmpparray,(np0+100))
      np0=np0+100  
    end if
    tmpparray(k)=present%prtcl   
    present=>present%next
  end do
  np=k
  parray=>allocateParticleArray(np)
  parray(1:np)=tmpparray(1:np)
  if (associated(tmpparray)) deallocate(tmpparray)
end subroutine toArray

end module linkedlist
