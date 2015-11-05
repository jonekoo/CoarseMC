!> A module for reading and writing coordinates of the particles and
!! the dimensions of the simulation box in a file in a general way. The
!! routines in the module don't define how the particle is written or
!! read but defines where it is read from or written to. Basicly it
!! defines the structure of the coordinate file but leaves details of
!! writing and reading to the particle and simulation box modules.
module class_factory
  use particle, only: particledat, readparticle, writeparticle
  use class_poly_box, only: poly_box, pbox_read, pbox_write
  use utils
  implicit none
  
  character(len=14), parameter :: beginmark = '#configuration'
  character(len=18), parameter :: endmark = '#end_configuration'
  
  type factory
     private
     integer :: unit = 5 
  end type factory

contains

  !> Writes the coordinates of @p particles and the dimensions of the 
  !! cylindrical simulation cell to the output unit specified by @p
  !! factory.
  !! 
  !! @param[in] afactory the object conducting the writing. 
  !! @param[in] outunit the Fortran io unit to write to. 
  !! @param[in] simbox the simulation box to be written.
  !! @param[in] particles the array of particles to write to the file.
  !! 
  subroutine factory_writestate(afactory, outunit, simbox, particles)
    type(factory), intent(in) :: afactory
    integer, intent(in) :: outunit
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    integer :: i
    !! Write simulation box.
    write(outunit, '(A)') beginmark
    call pbox_write(outunit, simbox)
    !! Write size of particle array to assist reading.
    write(outunit, '(/,' // fmt_char_int() //')') size(particles)
    !! Write particles
    do i = 1, size(particles) 
       call writeparticle(outunit, particles(i))
       write(outunit, '(A)') ''
    end do
    write(outunit, '(A)') endmark
  end subroutine factory_writestate
  
  !> Conducts the reading of the simulation box and particles from a file. 
  !! 
  !! @param[in] afactory the factory object that defines the order of reads.
  !! @param[in] inunit the Fortran io unit to read data from.
  !! @param[in,out] boxread the simulation box read. 
  !! @param[in,out] particles the particles read.
  !! @param[out] ios the status of the read operation at return. Negative value
  !! means an end of file condition and a positive value indicates an error as
  !! dictated in the Fortran 90/95 standard.
  !!
  subroutine factory_readstate(afactory, inunit, boxread, particles, ios)
    type(factory), intent(in) :: afactory
    integer, intent(in) :: inunit
    type(poly_box), intent(inout) :: boxread
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
    call pbox_read(inunit, boxread, ios)
    if (0 /= ios) return
    !! read dimension of particlearray
    read(inunit, *, iostat = ios) string
    read(string, *, iostat = ios) nparticles
    if (ios /= 0) then 
       write(*, *) 'Failed reading nparticles. Got ', string
       return
    end if
    
    if (allocated(particles)) then
       if (size(particles) /= nparticles) then
          deallocate(particles)
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
       call readparticle(inunit, particles(i), ios) 
       if (0 /= ios) return
    end do
    !! Read endmark
    read(inunit, *) string
    if (trim(string) /= endmark) then
       ios = 997
       return
    end if
  end subroutine factory_readstate

end module
