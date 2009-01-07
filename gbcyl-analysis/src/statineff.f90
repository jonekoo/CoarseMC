program statineff
implicit none
integer, parameter :: dp=selected_real_kind(12)
real(dp) :: blockvar
real(dp) :: totvar
real(dp) :: blockavg
real(dp) :: totavg
real(dp), dimension(:), allocatable :: values
real(dp) :: avalue, sumvalues
integer :: blocksize
integer :: astat, ios, i, valuecount
character(len=40) :: datafile, buf
integer, parameter :: maxsize=1000000
integer, parameter :: readunit=20
integer :: blockcount
  if (iargc()<2) stop 'Program needs arguments [file] [blocksize]'
  call getarg(1,buf)
  read(unit=buf,fmt=*) datafile 
  call getarg(2,buf)
  read(unit=buf,fmt=*) blocksize
  open(readunit, file=datafile, status='old', iostat=ios)
  if ( ios /= 0 ) then 
    write(*,*) 'Couldnt open file ',datafile
    stop;
  end if
  allocate(values(1:maxsize), stat=astat) 
  if(astat/=0) stop 'statineff: memory allocation failed'
  i=0
  do 
    read(unit=readunit,fmt=*,iostat=ios) avalue
    if (ios<0) exit;
    i=i+1
    if(i>=maxsize) stop 'the file is too big'
    values(i)=avalue
  end do  
  valuecount=i-mod(i,blocksize)
  !!sumvalues=sum(values(1:valuecount)); 
  !!totavg=sumvalues/valuecount
  !!totvar=0.0
  !!do i=1,valuecount
    !!totvar=totvar+(values(i)-totavg)**2
  !!end do
  !!totvar=totvar/valuecount
  !!blockvar=0.0
  !!blockcount=0
  !!do i=blocksize,valuecount,blocksize
    !!blockavg=sum(values(i-blocksize+1:i))/blocksize
    !!blockvar=blockvar+(blockavg-totavg)**2 
    !!blockcount=blockcount+1
  !!end do
  !!blockvar=blockvar/blockcount
  s_blockSize = blocksize*blockedVariance(values, valuecount, blocksize)/variance(values, valuecount) 
  write(*,*) s_blockSize
  contains 
    function blockedVariance(dataTable, nData, blockSize) result(blockedVar)
	  real::blockedVar=0.0
	  integer::
	end function blockedVariance
end program statineff


