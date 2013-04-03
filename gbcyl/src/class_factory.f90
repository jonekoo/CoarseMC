!> A module for reading and writing coordinates of particles and simulation 
!! box in a file in a general way. The routines in the module don't define how
!! the particle is written or read but defines where it is read from or written
!! to. Basicly it defines the structure of the coordinate file but leaves 
!! details of writing and reading to the particle and simulation box modules.
module class_factory
use particle
use class_poly_box
use utils
implicit none

character(len=14), parameter :: beginmark = '#configuration'
character(len=18), parameter :: endmark = '#end_configuration'

type factory
  private
  integer :: unit = 5 
end type

interface writestate
  module procedure factory_writestate
end interface

!interface binwritestate
!  module procedure factory_binwritestate
!end interface

interface readstate
  module procedure factory_readstate
end interface

contains

!> Writes the coordinates of @param particles and the dimensions of the 
!! cylindrical simulation cell to the file specified by @param factory.
!! 
!! @param factory the object conducting the writing. 
!! @param outunit the Fortran io unit to write to. 
!! @param simbox the simulation box to be written.
!! @param particles the array of particles to write to the file.
!! 
subroutine factory_writestate(afactory, outunit, simbox, particles)
  type(factory), intent(in) :: afactory
  integer, intent(in) :: outunit
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer :: i
  !! Write simulation box.
  write(outunit, '(A)') beginmark
  call write(outunit, simbox)
  !! Write size of particle array to assist reading.
  write(outunit, '(/,' // fmt_char_int() //')') size(particles)
  !! Write particles
  do i = 1, size(particles) 
    call write(outunit, particles(i))
    write(outunit, '(A)') ''
  end do
  write(outunit, '(A)') endmark
end subroutine

!> Writes the coordinates of @param particles and the dimensions of the 
!! cylindrical simulation cell to the file specified by @param factory.
!! 
!! @param factory the object conducting the writing. 
!! @param outunit the Fortran io unit to write to. 
!! @param simbox the simulation box to be written.
!! @param particles the array of particles to write to the file.
!! 
!subroutine factory_binwritestate(afactory, outunit, simbox, particles)
!  type(factory), intent(in) :: afactory
!  integer, intent(in) :: outunit
!  type(poly_box), intent(in) :: simbox
!  type(particledat), dimension(:), intent(in) :: particles
!  integer :: i
!  !! Write simulation box.
!  write(outunit) beginmark
!  call binwrite(outunit, simbox)
!  !! Write size of particle array to assist reading.
!  write(outunit) size(particles)
!  !! Write particles
!  do i = 1, size(particles) 
!    call binwrite(outunit, particles(i))
!    write(outunit) ''
!  end do
!  write(outunit) endmark
!end subroutine

!> Conducts the reading of the simulation box and particles from a file. 
!! 
!! @param afactory the factory object that defines the order of reads.
!! @param inunit the Fortran io unit to read data from.
!! @param boxread the simulation box read. 
!! @param particles the particles read.
!! @param iostatus the status of the read operation at return. Negative value
!! means an end of file condition and a positive value indicates an error as
!! dictated in the Fortran 90/95 standard.
!!
!! :TODO: Add iostatus variables to component reading routines such as the one 
!! :TODO: used to read individual particles?
!!
subroutine factory_readstate(afactory, inunit, boxread, particles, ios)
  type(factory), intent(in) :: afactory
  integer, intent(in) :: inunit
  type(poly_box), intent(out) :: boxread
  type(particledat), dimension(:), allocatable, intent(inout) :: particles
  integer, intent(out) :: ios
  integer :: nparticles
  integer :: astat
  integer :: i
  character(len = 80) :: string

  !! read beginmark
  read(inunit, *, iostat=ios) string
  if (ios /= 0) then
    return
  end if
  !! test validity of beginmark
  if (trim(string) /= beginmark) then
    ios = 998
    return
  end if
  !! read boxread
  call read(inunit, boxread, ios)
  if (0 /= ios) return
  !! read dimension of particlearray
  read(inunit, *, iostat = ios) string
  read(string, *, iostat = ios) nparticles
  if (ios /= 0) then 
    write(*, *) 'Failed reading nparticles. Got ', string
    return
  end if

  !! allocate space for particles.
  !if(associated(particles)) then
  if (allocated(particles)) then
    if (size(particles) /= nparticles) then
      deallocate(particles)
      !particles=>NULL()
    end if
  end if
  if( .not. allocated(particles)) then
    allocate(particles(nparticles), stat = astat)
    if (astat /= 0) then
      stop 'factory_readstate: Memory allocation for particles failed.'
    end if
  end if   
  
  !! read particles
  do i = 1, nparticles
    call read(inunit, particles(i), ios) 
    if (0 /= ios) return
  end do
  !! Read endmark
  read(inunit, *) string
  if (trim(string) /= endmark) then
    ios = 997
    return
  end if
end subroutine

end module
