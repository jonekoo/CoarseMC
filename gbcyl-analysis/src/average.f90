program average
implicit none
character(len=40) :: buf
character(len=40) :: inputfile
integer :: ios,nnumbers
integer, parameter :: input=20
real(kind=selected_real_kind(12)) :: anumber,sumnumbers
  if (iargc()<1) stop 'Program needs argument [file]'
  call getarg(1,buf)
  read(unit=buf,fmt=*) inputfile 
  open(input,file=inputfile,status='old',iostat=ios)
  if ( ios /= 0 ) then 
    write(*,*) 'Couldnt open file',inputfile
    stop;
  end if
  anumber=0.0
  sumnumbers=0.0
  nnumbers=0
  do 
    read(unit=input,fmt=*,iostat=ios) anumber
    if (ios<0) exit;
    sumnumbers=sumnumbers+anumber
    nnumbers=nnumbers+1
  end do  
  write(*,*) sumnumbers/real(nnumbers)
end program average
