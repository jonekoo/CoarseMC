!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
!! :TODO: make maximum number of bins in histogram and the interval between
!! bins adjustable!
!!
program pair_distribution
  use num_kind, only: dp
  use m_particledat, only : particledat
  use m_particle_factory, only: factory_readstate, factory
  use class_poly_box, only: poly_box
  use distribution
  use utils, only: fmt_char_dp
  use orientational_ordering
  implicit none
  type(particledat), allocatable :: particles(:)
  integer :: io_status
  integer, parameter :: maxbin = 100
  real(dp) :: delr = 0.1_dp
  real(dp), dimension(maxbin) :: histogram
  type(factory) :: afactory
  integer, parameter :: stdin = 5
  type(poly_box) :: simbox
  character(len=32) :: arg
  logical :: is_ljgb = .false., is_1D = .false., is_director = .false.
  logical, allocatable :: mask_i(:)
  integer :: i
  real(dp) :: p2
  logical :: is_r_bound = .false.

  do i = 1, command_argument_count()
     call get_command_argument(i, arg)

     select case (arg)
     case ('-z')
       !! 1D distributions are calculated with respect to z-axis of the
       !! simulation box.
       is_1D = .true.
     case ('-d')
       !! 1D distributions are calculated with respect to director.
       is_1D = .true.
       is_director = .true.
     case ('-h', '--help')
        call print_help()
        stop
     case ('--ljgb')
        is_ljgb = .true.
     case ('-r')
        is_r_bound = .true. 
     case default
        print '(a,a,/)', 'Unrecognized command-line option: ', arg
        call print_help()
        stop
     end select
  end do

  do  
    call factory_readstate(afactory, stdin, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    else
      allocate(mask_i(size(particles)))
      if (is_ljgb) then 
        mask_i = .not. particles%rod
      else
        mask_i = particles%rod
      end if
      if (is_1D) then
        direction = (/0._dp, 0._dp, 1._dp/)
        if (is_director) call orientation_parameter(pack(particles, &
          particles%rod), count(particles%rod), p2, direction)
        if (is_r_bound) then
          cutoff = 2._dp
          slice_area = 4._dp * atan(1._dp) * cutoff**2
          call distribution_func(simbox, particles, mask_i, particles%rod,&
          maxbin, delr, histogram, distance_1d_cyl, slice_volume)
        else 
          slice_area = simbox%volume() / getz(simbox)
          call distribution_func(simbox, particles, mask_i, particles%rod,&
          maxbin, delr, histogram, distance_1d, slice_volume)
        end if
      else
        call distribution_func(simbox, particles, mask_i, particles%rod,&
          maxbin, delr, histogram, distance_3d, spherical_shell_volume)
      end if
      write(*, '(100(' // fmt_char_dp() //',1X))') histogram
      if (allocated(mask_i)) deallocate(mask_i)
    end if
  end do

contains 

subroutine print_help()
  write(*, *) "USAGE:"
  write(*, *) "pair_distribution [OPTIONS] < molecule_file"
  write(*, *) ""
  write(*, *) "Without further options, pair_distribution prints the 3D radial"
  write(*, *) "pair distribution function of the Gay-Berne particles for each"
  write(*, *) "snapshots in the molecule_file"
  write(*, *) ""
  write(*, *) "OPTIONS:"
  write(*, *) ""
  write(*, *) "-z calculate 1D pair distribution in z-direction instead."
  write(*, *) "-d calculate 1D pair distribution with respect to the director instead."
  write(*, *) "-h, --help prints this help message."
  write(*, *) "--ljgb calculate pair distribution only for LJ-GB particle pairs." 
  write(*, *) ""
end subroutine

end program


