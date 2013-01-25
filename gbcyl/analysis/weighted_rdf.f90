!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
!! :TODO: make maximum number of bins in histogram and the interval between
!! bins adjustable!
!!
program analysis
  use nrtype, only: dp
  use particle, only : particledat
  use class_factory
  use class_poly_box
  use utils
  use gr3dweighted
  use layernormal
  implicit none
  type(particledat), dimension(:), pointer :: particles
  integer :: io_status
  integer, parameter :: maxbin = 100
  real(dp) :: delr = 0.1_dp
  real(dp), dimension(maxbin) :: weightedhistogram
  real(dp), dimension(maxbin) :: histogram
  type(factory) :: afactory
  integer, parameter :: stdin = 5
  type(poly_box) :: simbox
  real(dp), dimension(:, :), allocatable :: orientations
  integer :: i
  real(dp), parameter :: cutoff = 2._dp !! Change to non-hardcoded
  do  
    call readstate(afactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else
      if (.not. allocated(orientations)) then
        allocate(orientations(3, size(particles)))
      end if
      !write(*, *) size(particles)
      do i = 1, size(particles)
      !  write(*, *) i, orientation(particles(i)) 
        orientations(1:3, i) = localnormal(simbox, particles, i, cutoff) !orientation(particles(i))
      end do
      !forall(i = 1:size(particles))
      !  write(*, *) i, orientation(particles(i)) 
      !  orientations(1:3, i) = orientation(particles(i))
      !end forall
      call gr3dweighted1(particles, simbox, maxbin, delr, orientations, histogram, weightedhistogram)  
      write(*, '(100('//fmt_char_dp()//',1X))') weightedhistogram
    end if
  end do
  if(associated(particles)) deallocate(particles)
  !deallocate(orientations)
end program


