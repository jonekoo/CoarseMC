module utils
  use iso_fortran_env, only: iostat_end, iostat_eor
  use num_kind
  use mt_stream, only: rng=>genrand_double1_s, rngstate=>mt_state
  implicit none
  intrinsic index

  interface
     !> Computes the eigenvectors and eigenvalues @p values of @p matrix.
     !! @p matrix is replaced by the eigenvectors on output.
     subroutine eigensystem(matrix, values)
       use num_kind, only: dp
       implicit none
       real(dp), intent(inout) :: matrix(3, 3)
       real(dp), intent(out) :: values(3)
     end subroutine eigensystem
  end interface

  type str_wrapper
     character(len=:), allocatable :: c
  end type str_wrapper
    
contains 

!> Implements the Metropolis acceptance rule for a Monte Carlo update. 
!!
!! @param oldenergy the energy/enthalpy of the system before the move.
!! @param newenergy the energy/enthalpy of the system after the move.
!! @param genstate is the random number generator state.
!! @param isaccepted == .true. if the move is accepted. 
!!
  pure subroutine acceptchange(oldenergy, newenergy, temperature, genstate, &
       isaccepted)
  real(dp), intent(in) :: oldenergy, newenergy, temperature
  type(rngstate), intent(inout) :: genstate
  logical, intent(out) :: isaccepted
  real(dp) :: dE
  real(dp) :: r
  isaccepted = .true.
  dE = newenergy - oldenergy
  if (dE > 0._dp) then
     call rng(genstate, r)
     isaccepted = (r < exp(-dE/temperature))
  end if
end subroutine acceptchange

!> Returns the number of @p substr occurences in string @p str.
function substrcount(str, substr)
  character(len=*), intent(in) :: str
  character(len=*), intent(in) :: substr
  integer :: substrcount
  integer :: pos
  integer :: step
  substrcount = 0
  if (len(substr) == 0) then
    return
  else 
    pos = 1
    do while(.true.)
      step = index(str(pos:), substr)     
      if (step == 0) exit
      substrcount = substrcount + 1 
      pos = pos + step - 1 + len(substr)
    end do
  end if
end function


!> Returns a string @p joined that has been created by joining the
!! strings in @p str_array. Whitespace is trimmed away from elements of
!! str_array but not from @p separator.
subroutine join(str_array, separator, joined)
  character(len=*), intent(in) :: str_array(:)
  character(len=*), intent(in) :: separator
  character(len=*), intent(inout) :: joined
  integer :: i, joined_len
  joined_len = 0
  !! Calculate the length of the joined string
  do i = 1, size(str_array)
    joined_len = joined_len + len_trim(str_array(i))
  end do
  joined_len = joined_len + (size(str_array) - 1) * len(separator)
  if (len(joined) < joined_len) stop 'utils: join: Error given string joined is too short.' 
  joined = trim(adjustl(str_array(1))) 
  do i = 2, size(str_array)
    joined = trim(adjustl(joined)) // separator // trim(adjustl(str_array(i)))
  end do  
end subroutine


!> Splits a given string @p str to the array @p strarr using @p 
!! delimiterstr as the delimiter. Delimiters are not included in the results.
!! Result may contain empty strings. 
subroutine splitstr(str, delimiterstr, strarr)
  character(len=*), intent(in) :: str
  character(len=*), intent(in) :: delimiterstr
  character(len=*), dimension(:), allocatable, intent(inout) :: strarr
  integer :: dim, i 
  integer :: lastpos
  integer :: step
  dim = substrcount(str, delimiterstr) + 1
  allocate(strarr(dim))
  lastpos = 1
  do i = 1, dim
    step = index(str(lastpos:), delimiterstr)
    if (step /= 0) then
      strarr(i) = str(lastpos:lastpos + step - 2)
    else ! last substr
      strarr(i) = str(lastpos:)
    end if
    lastpos = lastpos + step - 1 + len(delimiterstr)
  end do
end subroutine

!> Returns the components of the vector @p orientation in cylindrical
!! coordinates when the vector is at @p position.
pure function unitvec(orientation, position)
  real(dp), dimension(3), intent(in) :: orientation
  real(dp), dimension(3), intent(in), optional :: position
  real(dp), dimension(3) :: unitvec
  real(dp) :: uro, utheta, uz
  real(dp) :: nx, ny, nz, theta
  theta = atan2(position(2), position(1))
  nx = 0._dp
  ny = 0._dp
  nz = 1._dp
  call rotate_vector(orientation(1), orientation(2), orientation(3), nx, ny, nz, &
  theta, uro, utheta, uz)
  unitvec = (/uro, utheta, uz/)
end function
  

!> Rotates the vector (@p x, @p y, @p z) into (@p xp, @p yp, @p zp)
!! around axis (@p nx, @p ny, @p nz) [unit vector of the direction]
!! through angle @p angle. 
!!
!! Author Juho Lintuvuori.
!! Modified by Jouni Karjalainen to kind parameter dp and intrinsic dot_product
!!
!! The angle is defined counterclockwise when looking at the direction of 
!! (nx, ny, nz).
!!
!! @see Goldstein: Classical Mechanics 2nd ed., p. 165
!!
pure subroutine rotate_vector(x, y, z, nx, ny, nz, angle, xp, yp, zp)
  implicit none
  real(dp), intent(in) :: x, y, z, nx, ny, nz, angle
  real(dp), intent(out) :: xp, yp, zp
  real(dp) :: dotpr, cp, sp
  real(dp), dimension(3) :: vec1 
  real(dp), dimension(3) :: vec2 
  vec1 = (/nx, ny, nz/)
  vec2 = (/x, y, z/)
  dotpr = dot_product(vec1, vec2)
  cp = cos(angle)
  sp = sin(angle)
  call crossp(x, y, z, nx, ny, nz, xp, yp, zp)
  xp = x * cp + nx * dotpr * (1._dp - cp) + xp * sp
  yp = y * cp + ny * dotpr * (1._dp - cp) + yp * sp
  zp = z * cp + nz * dotpr * (1._dp - cp) + zp * sp
end subroutine rotate_vector



!> Rotates a rank 2 tensor defined in @p old_axes to @p new_axes. 
!!
!! @param original the tensor to be rotated.
!! @param old_axes the axes in which the original tensor is defined. 
!! @param new_axes the axes in which we want to represent the tensor.
!! @param rotated the @p original tensor rotated from @p old_axes to @p new_axes.
!!
subroutine rotate_tensor(original, old_axes, new_axes, rotated)
  real(dp), intent(in) :: original(3, 3), old_axes(3, 3), new_axes(3, 3)
  real(dp), intent(out) :: rotated(3, 3)

  real(dp) :: cosines(3, 3)
  integer :: a, b, i, j

  do a = 1, 3
    do i = 1, 3
      cosines(a, i) = dot_product(new_axes(1:3, a), old_axes(1:3, i)) 
    end do
  end do

  rotated = 0._dp
  do a = 1, 3
  do b = 1, 3
    do i = 1, 3
      do j = 1, 3
        rotated(a, b) = rotated(a, b) + cosines(a, i) * cosines(b, j) * original(i, j)
      end do
    end do
  end do
  end do
end subroutine  



!> Returns the cross-product c = a x b.
pure function cross_product(a, b) result(c)
  real(dp), intent(in) :: a(3), b(3)
  real(dp) :: c(3)
  c = 0._dp
  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = -a(1) * b(3) + a(3) * b(1)
  c(3) = a(1) * b(2) - a(2) * b(1)
end function

!> Computes the cross-product (@p CX, @p CY, @p CZ) =
!! (@p AX, @p AY, @p AZ) x (@p BX, @p BY, @p BZ) 
pure SUBROUTINE CROSSP(AX,AY,AZ,BX,BY,BZ,CX,CY,CZ)
  REAL(DP), intent(in) :: AX, AY, AZ, BX, BY, BZ
  real(dp), intent(out) :: CX, CY, CZ
  CX = AY * BZ - AZ * BY
  CY = AZ * BX - AX * BZ
  CZ = AX * BY - AY * BX
END SUBROUTINE

!> Returns the value of a polynomial of order n-1 in the point @p x. 
!! Coefficients of the polynomial are given in table @p a, where 
!! @p a(i) is the coefficient of the x^(n-i) term. Calculation is done
!! with the Horner method. Size of @p a is n.
pure function horner(a, x) 
  real(dp), dimension(:), intent(in) :: a
  real(dp), intent(in) :: x 
  real(dp) :: horner
  integer :: i
  horner = a(1)
  do i = 2, size(a)
    horner = horner * x + a(i)
  end do 
end function horner


!> Returns a formatting character for the default integer type. To be
!! used when consistent formatting for integer output is needed.
pure function fmt_char_int() result(format_char)
  character(len = 50) :: format_char
  integer :: r
  character(len = 50) :: w_char
  r = range(r)
  write(w_char, *) 1 + r !! sign + range
  write(format_char, *) 'I' // trim(adjustl(w_char))
  format_char = trim(adjustl(format_char))
end function


!> Returns a formatting character for a double precision real number.
!! To be used when consistent formatting of real numbers is needed. 
pure function fmt_char_dp() result(format_char)
  integer :: e
  integer :: w
  integer :: d
  real(dp) :: u 
  character(len = 50) :: w_char !! Width of field
  character(len = 50) :: d_char !! Width of decimal field
  character(len = 50) :: e_char !! width of exponent
  character(len = 50) :: format_char
  e = int(log10(real(range(u)))) + 1
  d = precision(u)
  w = 3 + d + 2 + e
  write(w_char, *) w
  write(e_char, *) e
  write(d_char, *) d
  write(format_char, *) 'G' // trim(adjustl(w_char)) // '.' // &
       trim(adjustl(d_char)) // 'E' // trim(adjustl(e_char)) 
  format_char = trim(adjustl(format_char))
end function

!> Returns a formatting character for an array of double precision real
!! numbers. To be used when consistent formatting of real numbers is
!! needed. 
pure function fmt_char_dp_array(array_size) result(format_char)
  integer, intent(in) :: array_size
  character(len = 50) :: format_char
  write(format_char, *) '(', array_size - 1, '(' // &
       trim(adjustl(fmt_char_dp())) // ',1X),' // &
       trim(adjustl(fmt_char_dp())) // ')'
  format_char = trim(adjustl(format_char))
end function fmt_char_dp_array

subroutine readstr(unit, str, ios)
  integer, intent(in) :: unit
  character(len=:), allocatable, intent(out) :: str
  integer, intent(out) :: ios
  character(len=1) :: buf
  str = ''
  do while(.true.)
     read(unit=unit, advance='no', iostat=ios, fmt='(A)') buf
     !write(*, *) buf
     if (ios == iostat_end) then
        return
     else if (ios == iostat_eor) then
        str = str // new_line(str)
        continue
     else if (ios /= 0) then
        return
     else
        str = str // buf
     end if
  end do
end subroutine readstr
  
subroutine readlines(unit, lines, ios)
  intrinsic new_line
  integer, intent(in) :: unit
  type(str_wrapper), allocatable, intent(out) :: lines(:)
  integer, intent(out) :: ios
  character(len=:), allocatable :: temp
  integer :: current, next
  call readstr(unit, temp, ios)
  allocate(lines(0))
  current = 1
  do while(current < len(temp))
     next = index(temp(current:), new_line(temp))
     if (next == 0) then
        lines = [lines, str_wrapper(temp(current:))]
        exit
     else
        lines = [lines, str_wrapper(temp(current:next))]
     end if
     current = next + 1
  end do
end subroutine readlines

end module utils
