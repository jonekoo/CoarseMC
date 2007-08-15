  !Kirjoittaa tilan eli sylinterin säteen, korkeuden ja partikkelitaulukon
  !tiedostoon  
  subroutine writestate(T,P,R,Lz,particlearray,index)
    implicit none
    intrinsic trim,adjustl
    real(dp), intent(in) :: T,P,R,Lz
    integer, intent(in) :: index
    type(particledat),dimension(:),pointer :: particlearray
    character(len=30) :: outfile
    integer :: GB=0,Xe=0,ios,N,astat,i
    integer, parameter :: outunit=19
    integer, dimension(:), allocatable :: help
    N=size(particlearray)
    allocate(help(N),stat=astat)
    if (astat/=0) then 
      write(*,*) 'writestate: Virhe varattaessa muistia taulukolle help.'
      stop
    end if
    GB=0
    Xe=0
    do i=1,N
      if(particlearray(i)%rod) then 
        GB=GB+1
        help(i)=1
      else 
        Xe=Xe+1
        help(i)=0
      end if
    end do 
    outfile=trim(adjustl(filename(R,T,P,GB,Xe,index)) )
    outfile=trim(outfile)
    open(outunit,file=outfile,status='replace',position='append',form='formatted',iostat=ios)
    if (ios/=0) then
      write (*,*) 'writestate: tiedostoa ',outfile,'ei voitu avata.'
      stop;
    end if
    write(outunit,*) '$R:',R,'$Lz:',Lz
    write(outunit,*) '$N:',N,'$GB:',GB,'$Xe:',Xe
    write(outunit,*) '$x:'
    write(outunit,*) particlearray%x
    write(outunit,*) '$y:'
    write(outunit,*) particlearray%y
    write(outunit,*) '$z:'
    write(outunit,*) particlearray%z
    write(outunit,*) '$rod:'
    write(outunit,*) help
    write(outunit,*) '$ux:'
    write(outunit,*) particlearray%ux
    write(outunit,*) '$uy:'
    write(outunit,*) particlearray%uy
    write(outunit,*) '$uz:'
    write(outunit,*) particlearray%uz
    close(outunit) 
    deallocate(help)
  end subroutine writestate
