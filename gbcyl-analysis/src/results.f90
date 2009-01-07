module results

  character(len = *), parameter :: filename_orientation = 'eigens.dat'   
  integer, parameter :: unit_orientation = 10  
  logical, save :: is_initialized = .false.
  integer, parameter :: unit_tau1 = 11
  character(len = *), parameter :: filename_tau1 = 'tau1.dat'



  contains 



  subroutine init
    implicit none
    call open_file(filename_orientation, unit_orientation)
    call open_file(filename_tau1, unit_tau1)
    !! Write labels to resultFiles
    is_initialized = .true. 
  end subroutine init



  subroutine finalize
    implicit none
    close(unit_orientation)
    close(unit_tau1)
    is_initialized = .false.
  end subroutine finalize



  subroutine open_file(fileName, writeUnit)
    implicit none
    character(len = *), intent(in) :: fileName
    integer, intent(in) :: writeUnit
    integer :: fileStatus
    open(writeUnit, FILE = fileName, STATUS = 'NEW',&
       & POSITION = 'APPEND', IOSTAT = fileStatus)
    if(fileStatus > 0) then
      write(*, *) 'Error ', fileStatus, 'when opening file ', fileName 
      write(*, *) 'Program will end.'
      stop
    end if
  end subroutine open_file



  subroutine calculate(particleArray, n_particles, radius, height)
    use particle, only: particledat
    use nrtype
    use orientational_ordering, only: eigens
    use translational_ordering, only: tau1
    implicit none
    intrinsic MAXLOC
    
    type(particledat), dimension(:), pointer :: particleArray
    integer, intent(in) :: n_particles
    real(dp), intent(in) :: radius, height
    integer, parameter :: n_dimensions = 3
    real(dp), dimension(n_dimensions) :: values
    real(dp), dimension(n_dimensions, n_dimensions) :: vectors
    character(len = 500) :: entry_orientation
    integer, dimension(1) :: max_value_pos
    real(dp), dimension(3) :: direction
    real(dp) :: tau1_value
    real(dp) :: layer_distance
    character(len = 500) :: entry_tau1

    if (.not. is_initialized) then
      write(*, *) 'Tried to use calculate before initializing module.'
      write(*, *) 'Program will end.'
      stop
    end if

    !! Calculate eigenvalues and eigenvectors of the orientational
    !! ordering tensor. 
    call eigens(particleArray, n_particles, values, vectors)
    write(entry_orientation, *) values, vectors
    write(unit_orientation, *) TRIM(entry_orientation)

    !! Get the eigenvector corresponding to the largest eigenvalue
    !! and calculate the translational order parameter using that 
    !! direction
    max_value_pos = MAXLOC(values)
    direction(1:n_dimensions) = vectors(max_value_pos(1), 1:n_dimensions)
    call tau1(particleArray, n_particles, direction, tau1_value, layer_distance)
    write(entry_tau1, *) tau1_value, layer_distance, direction
    write(unit_tau1, *) TRIM(entry_tau1)

  end subroutine calculate



end module results
