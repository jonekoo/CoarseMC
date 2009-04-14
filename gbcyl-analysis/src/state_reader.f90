module state_reader
  use nrtype, only: dp
  use particle, only: particledat


   
  contains



  subroutine open_read_unit(fileName, read_unit)
    implicit none
    character(len = *), intent(in) :: filename
    integer, intent(in) :: read_unit
    integer :: ios
    open(read_unit, file = filename, status = 'old', iostat = ios)
    if (ios /= 0) then
      stop 'Opening of file failed.'
    end if
  end subroutine open_read_unit


!!  subroutine read_configuration(particles, n_particles, radius, height, & 
!!    iostatus)
!  implicit none
!  type(particledat), dimension(:), pointer :: particles
!  integer, intent(out) :: n_particles
!  real(dp), intent(out) :: radius, height
!  integer, intent(out) :: iostatus
!  integer :: read_unit = 5
!    call read_configuration(particles, n_particles, radus, height, iostatus, &
!      read_unit)
!  end subroutine read_configuration


  subroutine read_configuration(particles, n_particles, radius, height, &
    iostatus, runit)
  implicit none
  type(particledat), dimension(:), pointer :: particles
  integer, intent(out) :: n_particles
  real(dp), intent(out) :: radius, height
  integer, intent(out) :: iostatus
  integer, intent(in), optional :: runit 
    integer :: particlesstat, helpStat, i
    character(len = 3) :: charvar
    integer, dimension(:), allocatable :: help
    logical :: readOpened
    integer :: read_unit
    if(PRESENT(runit)) then
      read_unit = runit
    else 
      read_unit = 5
    end if
    radius = 0.0
    height = 0.0 
    iostatus = 0
    inquire(UNIT = read_unit, OPENED = readOpened)
    if (.not. readOpened) stop 'read_unit not open'
    read(read_unit, *, IOSTAT = iostatus) charvar, radius, charvar, height
    read(read_unit, *, IOSTAT = iostatus) charvar, n_particles
    read(read_unit, *, IOSTAT = iostatus) charvar
    particlesstat = 0
    if (size(particles) /= n_particles) then
      if (associated(particles)) deallocate(particles)
      allocate(particles(n_particles), stat = particlesstat)
    end if
    helpStat = 0
    allocate(help(n_particles), stat = helpStat)
    if (particlesstat /= 0 .or. helpStat /= 0 ) then
      write(*, *) 'readstate: particles, help'
      stop;
    end if
    read(read_unit, *, IOSTAT = iostatus) particles(1:n_particles)%x
    read(read_unit, *, IOSTAT = iostatus) charvar
    read(read_unit, *, IOSTAT = iostatus) particles(1:n_particles)%y
    read(read_unit, *, IOSTAT = iostatus) charvar, &
      particles(1:n_particles)%z
    read(read_unit, *, IOSTAT = iostatus) charvar
    read(read_unit, *, IOSTAT = iostatus) help(1:n_particles)
    read(read_unit, *, IOSTAT = iostatus) charvar
    read(read_unit, *, IOSTAT = iostatus) particles(1:n_particles)%ux
    read(read_unit, *, IOSTAT = iostatus) charvar, &
      particles(1:n_particles)%uy
    read(read_unit, *, IOSTAT = iostatus) charvar, &
      particles(1:n_particles)%uz
    do i=1,n_particles
      if(help(i)==1) then
        particles(i)%rod=.true.
      else
        particles(i)%rod=.false.
      end if
    end do
    deallocate(help)
  end subroutine read_configuration



  subroutine readerror(filename)
    implicit none
    character(len=*) :: filename
    write(*,*) 'Tiedoston ',filename,' muoto on v‰‰r‰.'
    write(*,*) 'Tilaa ei voitu lukea. Ohjelman suoritus keskeytyy'
    stop;
  end subroutine readerror

end module state_reader
