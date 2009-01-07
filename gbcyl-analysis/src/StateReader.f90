module StateReader
  use nrtype
  use particle

  integer, parameter :: readunit=17   
   
  contains

  subroutine openFile(fileName)
    implicit none
    character(len = *), intent(in) :: filename
    integer :: ios
    !! Koittaa avata annetun tiedoston
    open(readunit, file = filename, status = 'old', iostat = ios)
    if (ios /= 0) then
      write(*, *) 'Tiedoston ',filename,' avaaminen ep‰onnistui.'
      write(*, *) 'Ohjelman suoritus keskeytyy.'
      stop;
    end if
  end subroutine openFile



  subroutine closeFile
    implicit none
    close(readunit)
  end subroutine closeFile



  subroutine readConfiguration(particleArray, nParticles, R, Lz, iostatus)
    implicit none

    type(particledat), dimension(:), pointer :: particleArray
    real(dp), intent(out) :: R, Lz
    integer, intent(out) :: iostatus
    integer, intent(out) :: nParticles
    integer :: ios, particleArraystat, helpStat, i
    character(len = 3) :: charvar
    integer, dimension(:), allocatable :: help
    logical :: readOpened

    R = 0.0
    Lz = 0.0 
    iostatus = 0
    inquire(UNIT = readunit, OPENED = readOpened)
    if (.not. readOpened) stop 'readunit not open'

    read(readunit, *, IOSTAT = iostatus) charvar, R, charvar, Lz
    read(readunit, *, IOSTAT = iostatus) charvar, nParticles
    read(readunit, *, IOSTAT = iostatus) charvar
    !! Muistin varaus partikkelitaulukolle
    particleArraystat = 0
    if (size(particleArray) /= nParticles) then
      if (associated(particleArray)) deallocate(particleArray)
      allocate(particleArray(nParticles), stat = particleArraystat)
    end if
    helpStat = 0
    allocate(help(nParticles), stat = helpStat)
    if (particleArraystat /= 0 .or. helpStat /= 0 ) then
      write(*, *) 'readstate:Virhe varattaessa muistia: particleArray, help'
      stop;
    end if
      
    read(readunit, *, IOSTAT = iostatus) particleArray(1:nParticles)%x
    read(readunit, *, IOSTAT = iostatus) charvar
    read(readunit, *, IOSTAT = iostatus) particleArray(1:nParticles)%y
    read(readunit, *, IOSTAT = iostatus) charvar, particleArray(1:nParticles)%z
    read(readunit, *, IOSTAT = iostatus) charvar
    read(readunit, *, IOSTAT = iostatus) help(1:nParticles)
    read(readunit, *, IOSTAT = iostatus) charvar
    read(readunit, *, IOSTAT = iostatus) particleArray(1:nParticles)%ux
    read(readunit, *, IOSTAT = iostatus) charvar, particleArray(1:nParticles)%uy
    read(readunit, *, IOSTAT = iostatus) charvar, particleArray(1:nParticles)%uz
 
    do i=1,nParticles
      if(help(i)==1) then
        particleArray(i)%rod=.true.
      else
        particleArray(i)%rod=.false.
      end if
    end do
    deallocate(help)
  end subroutine readConfiguration



  subroutine readerror(filename)
    implicit none
    character(len=*) :: filename
    write(*,*) 'Tiedoston ',filename,' muoto on v‰‰r‰.'
    write(*,*) 'Tilaa ei voitu lukea. Ohjelman suoritus keskeytyy'
    stop;
  end subroutine readerror

end module StateReader
