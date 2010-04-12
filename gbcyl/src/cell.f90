module cell
  use nrtype
  implicit none
  private

  public :: list
  public :: iterator
  public :: new_list
  public :: ncells
  public :: cellindex
  public :: isdone
  public :: advance
  public :: value
  public :: new_iterator
  public :: nx, ny, nz
  public :: ix, iy, iz  
  public :: delete
  public :: writetostdout

  !! nx, ny and nz are the numbers of cells in x, y and z directions.
  !! ncells is the total number of cells 
  !! heads contains the linked list heads for each cell, 
  !! size(heads) == ncells. 
  !! list contains the one directional links between particles, 
  !! size(links) == npositions. 
  !! 
  type list
    private
    integer :: nx = 0 
    integer :: ny = 0
    integer :: nz = 0
    integer, dimension(:), pointer :: heads => NULL()
    integer, dimension(:), pointer :: links => NULL()
  end type list

  type iterator
    private
    type(list), pointer :: list => NULL()
    integer :: current = 0
  end type iterator

  interface delete
    module procedure cl_delete
  end interface

  interface new_list
    module procedure new_listf
  end interface

  interface ncells
     module procedure ncellsinlist, ncellsonside
  end interface

  interface cellindex
    module procedure indexr, indexi
  end interface
 
  interface isdone
    module procedure cit_isdone
  end interface
  
  interface advance
    module procedure cit_advance
  end interface
 
  interface value
    module procedure cit_value
  end interface

  interface nx
    module procedure cl_nx
  end interface

  interface ny
    module procedure cl_ny
  end interface

  interface nz
    module procedure cl_nz
  end interface

  interface assignment(=)
    module procedure assignlist
  end interface

contains

  pure subroutine assignlist(to, from)
    type(list), intent(in) :: from
    type(list), intent(inout) :: to
    integer :: nheads, nlinks
    if (associated(to%heads)) deallocate(to%heads)
    if (associated(to%links)) deallocate(to%links)
    to%nx = from%nx
    to%ny = from%ny
    to%nz = from%nz
    nheads = size(from%heads)
    nlinks = size(from%links)
    allocate(to%heads(size(from%heads)))
    to%heads(1:nheads) = from%heads(1:nheads) 
    allocate(to%links(size(from%links)))
    to%links(1:nlinks) = from%links(1:nlinks)
  end subroutine

  !! :TODO: add possibility to make an even number of cells in each direction
  !! :TODO: otherwise this function does not make sense.
  elemental function ncellsonside(boxside, minlength) result(n)
    real(dp), intent(in) :: boxside
    real(dp), intent(in) :: minlength
    integer :: n
    n = int(boxside / minlength)
  end function

  elemental function ncellsinlist(cl) result(n)
    type(list), intent(in) :: cl
    integer :: n
    n = cl%nx * cl%ny * cl%nz
  end function

  !! Constructs a new cell list. 
  !! 
  !! @see Allen, M. P. and Tildesley, D. J.: Computer simulation of liquids: 
  !! Neighbour lists.
  !! 
  !! @p positions are the positions of the objects added to the list scaled 
  !! so that the cartesian coordinate in each direction belongs to the 
  !! interval [-1, 1].
  !! @p nx number of cells in x-direction
  !! @p ny number of cells in y-direction
  !! @p nz number of cells in z-direction
  !!
  pure function new_listf(positions, nxin, nyin, nzin) result(cl)
    real(dp), dimension(:, :), intent(in) :: positions
    integer, intent(in) :: nxin, nyin, nzin
    type(list) :: cl
    integer :: ipos
    integer :: allocstat
    integer :: icell
    allocate(cl%links(size(positions)/3), stat = allocstat)
    !if (allocstat /= 0) then
    !  stop 'Memory allocation failed in new_listnew'
    !end if
    cl%nx = nxin 
    cl%ny = nyin 
    cl%nz = nzin 
    allocate(cl%heads(ncells(cl)), stat = allocstat)
    !if (allocstat /= 0) then
    !  stop 'Memory allocation failed in new_celllist'
    !end if
    do icell = 1, ncells(cl)
      cl%heads(icell) = 0
    end do
    do ipos = 1, size(positions)/3
      icell = cellindex(cl, positions(1:3, ipos))
      cl%links(ipos) = cl%heads(icell)
      cl%heads(icell) = ipos
    end do
  end function

  subroutine writetostdout(cl)
    type(list), intent(in) :: cl
    write(*, *) 'nx = ', cl%nx, 'ny = ', cl%ny, 'nz = ', cl%nz
    write(*, *) 'heads = ', cl%heads(1:cl%nx*cl%ny*cl%nz)
    write(*, *) 'links = ', cl%links
  end subroutine

  pure subroutine cl_delete(cl)
    type(list), intent(inout) :: cl
    if (associated(cl%heads)) deallocate(cl%heads)
    if (associated(cl%links)) deallocate(cl%links)
  end subroutine

  pure function indexr(cl, r) result(cellindex)
    type(list), intent(in) :: cl
    real(dp), dimension(3), intent(in) :: r
    integer :: cellindex
    cellindex = 1 + int((r(1) + 0.5_dp) * real(cl%nx, dp)) &
    + int((r(2) + 0.5_dp) * real(cl%ny, dp)) * cl%nx &
    + int((r(3) + 0.5_dp) * real(cl%nz, dp)) * cl%nx * cl%ny
  end function 

  elemental function ix(cl, x)
    type(list), intent(in) :: cl
    real(dp), intent(in) :: x
    integer :: ix
    ix = int((x + 0.5_dp) * real(cl%nx, dp))
  end function

  elemental function iy(cl, y)
    type(list), intent(in) :: cl
    real(dp), intent(in) :: y
    integer :: iy
    iy = int((y + 0.5_dp) * real(cl%ny, dp))
  end function

  elemental function iz(cl, z)
    type(list), intent(in) :: cl
    real(dp), intent(in) :: z
    integer :: iz
    iz = int((z + 0.5_dp) * real(cl%nz, dp))
  end function

  elemental function indexi(alist, ixpos, iypos, izpos) result(icell)
    implicit none
    integer :: icell
    type(list), intent(in) :: alist
    integer, intent(in) :: ixpos
    integer, intent(in) :: iypos
    integer, intent(in) :: izpos
    logical :: isxvalid
    logical :: isyvalid
    logical :: iszvalid
    isxvalid = ixpos >= 0 .and. ixpos < cl_nx(alist)
    isyvalid = iypos >= 0 .and. iypos < cl_ny(alist)
    iszvalid = izpos >= 0 .and. izpos < cl_nz(alist)
    if (isxvalid .and. isyvalid .and. iszvalid) then
      icell = 1 + mod(ixpos + alist%nx, alist%nx) &
      + mod(iypos  + alist%ny, alist%ny) * alist%nx &
      + mod(izpos  + alist%nz, alist%nz) * alist%nx * alist%ny 
    else 
      icell = 0
    end if
  end function

  elemental function cl_nx(alist)
    integer :: cl_nx
    type(list), intent(in) :: alist
    cl_nx = alist%nx
  end function

  elemental function cl_ny(alist)
    integer :: cl_ny
    type(list), intent(in) :: alist
    cl_ny = alist%ny
  end function

  elemental function cl_nz(alist)
    integer :: cl_nz
    type(list), intent(in) :: alist
    cl_nz = alist%nz
  end function
 
  !! Function which returns a mask corresponding to cell position 
  !! Particle indices in the cell will be marked .true.
  !!
  pure function maski(cl, position, n)
    type(list), intent(in) :: cl
    real(dp), dimension(3), intent(in) :: position
    integer, intent(in) :: n
    logical, dimension(n) :: maski
    integer :: index
    maski(:) = .false.
    index = cl%heads(cellindex(cl, position))
    do
      if (index == 0) exit
      maski(index) = .true.
      index = cl%links(index)
    end do  
  end function



  !! Iterator subroutines

  !! Returns a new iterator for the cell @p icell in @p alist.
  !!
  !! Would it be reasonable to write an iterator to go through all the cells?
  !! Maybe not since the list is probably in the scope of the user of the 
  !! iterator.
  !!
  function new_iterator(alist, icell)
    type(iterator) :: new_iterator
    type(list), intent(in), target :: alist
    integer, intent(in) :: icell
    if (icell < 1 .or. icell > ncells(alist)) then
      new_iterator%list => alist
      new_iterator%current = 0
    else
      new_iterator%list => alist
      new_iterator%current = alist%heads(icell)
    end if
  end function

  !! Advances the iterator in the cell list. 
  !! 
  !! There's no next(aniterator) function returning an iterator since it's
  !! not clear to the author how the list pointer assignment would work in 
  !! that case.
  !!
  elemental subroutine cit_advance(aniterator)
    type(iterator), intent(inout) :: aniterator
    if (.not. isdone(aniterator)) then
      aniterator%current = aniterator%list%links(aniterator%current)
    end if
  end subroutine

  elemental function cit_value(aniterator) result(current)
    integer :: current
    type(iterator), intent(in) :: aniterator
    current = aniterator%current
  end function

  elemental function cit_isdone(aniterator) result(done)
    logical :: done
    type(iterator), intent(in) :: aniterator
    done = (0 == aniterator%current)    
  end function

end module
