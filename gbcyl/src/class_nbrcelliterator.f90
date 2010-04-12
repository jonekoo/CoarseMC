module class_nbrcelliterator
implicit none
private
 
public :: value
public :: advance
public :: isdone
public :: new_nbrcelliterator1d
public :: new_nbrcelliterator
public :: nbrcelliterator
public :: nbrcelliterator1d
public :: nbrcellit_xvalue
public :: nbrcellit_yvalue
public :: nbrcellit_zvalue


interface value
  module procedure nbrcellit1d_value, nbrcellit_value
end interface

interface advance
  module procedure nbrcellit_advance, nbrcellit1d_advance
end interface

interface isdone
  module procedure nbrcellit1d_isdone, nbrcellit_isdone
end interface

type nbrcelliterator1d
  private
  integer :: n
  integer :: begin
  integer :: end
  integer :: current
  logical :: isdone
end type

type nbrcelliterator
  private
  integer :: ix, iy, iz
  integer :: nx, ny, nz
  type(nbrcelliterator1d) :: xit, yit, zit
end type

contains

pure function new_nbrcelliterator1d(n, i) result(it)
  integer, intent(in) :: n
  integer, intent(in) :: i
  type(nbrcelliterator1d) :: it
  if (i < 0 .or. i >= n) then
    it%isdone = .true.
  else
    it%isdone = .false.
  end if
  it%n = n
  it%begin = mod(i - 1 + it%n, it%n)
  it%end = mod(i + 1, it%n)
  if (n < 3) it%begin = i
  if (n < 2) it%end = i
  it%current = it%begin
end function

pure function nbrcellit1d_value(it) result(val)
  type(nbrcelliterator1d), intent(in) :: it
  integer :: val
  val = it%current
end function

pure function nbrcellit1d_isdone(it) result(done)
  type(nbrcelliterator1d), intent(in) :: it
  logical :: done
  done = it%isdone
end function

pure subroutine nbrcellit1d_advance(it)
  type(nbrcelliterator1d), intent(inout) :: it
  if (it%current == it%end) then
    it%isdone = .true.
  else 
    it%current = mod(it%current + 1, it%n)
  end if
end subroutine

pure function new_nbrcelliterator(nx, ny, nz, ix, iy, iz) result(it)
  integer, intent(in) :: nx, ny, nz
  integer, intent(in) :: ix, iy, iz
  type(nbrcelliterator) :: it
  it%ix = ix
  it%iy = iy
  it%iz = iz
  it%nx = nx
  it%ny = ny
  it%nz = nz
  it%xit = new_nbrcelliterator1d(nx, ix)
  it%yit = new_nbrcelliterator1d(ny, iy)
  it%zit = new_nbrcelliterator1d(nz, iz)
end function

pure function nbrcellit_isdone(it) result(done)
  type(nbrcelliterator), intent(in) :: it
  logical :: done
  done = isdone(it%zit) .or. isdone(it%yit) .or. isdone(it%xit)
end function

pure subroutine nbrcellit_advance(it)
  type(nbrcelliterator), intent(inout) :: it
  call advance(it%xit)
  if (isdone(it%xit)) then
    it%xit = new_nbrcelliterator1d(it%nx, it%ix)
    call advance(it%yit)
    if (isdone(it%yit)) then
      it%yit = new_nbrcelliterator1d(it%ny, it%iy)
      call advance(it%zit)
    end if
  end if
end subroutine

pure function nbrcellit_value(it) result(val)
  type(nbrcelliterator), intent(in) :: it
  integer :: val
  val = value(it%xit) + value(it%yit) * it%nx + value(it%zit) * it%nx * it%ny
end function

pure function nbrcellit_xvalue(it) result(xval)
  type(nbrcelliterator), intent(in) :: it
  integer :: xval 
  xval = value(it%xit)
end function

pure function nbrcellit_yvalue(it) result(yval)
  type(nbrcelliterator), intent(in) :: it
  integer :: yval 
  yval = value(it%yit)
end function

pure function nbrcellit_zvalue(it) result(zval)
  type(nbrcelliterator), intent(in) :: it
  integer :: zval 
  zval = value(it%zit)
end function

end module
