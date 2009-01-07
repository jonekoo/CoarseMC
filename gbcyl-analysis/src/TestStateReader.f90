

program TestStateReader
  implicit none
!! Tulostetaan toisen konfiguraation ensimmainen partikkeli
!! Luetaan konfiguraatiot tiedostosta ja tulostetaan niiden maara. 
 
  call testRods
  call testNumberOfConfigurations
contains 



  subroutine testRods
    use StateReader
    use particle
    use nrtype
    implicit none
    type(particledat), dimension(:), pointer :: particleArray
    real(dp) :: radius, height
    integer :: iostatus
    character(len = *), parameter :: fileName = 'simdata.out'
    integer :: i, j, n

    radius = 0.0
    height = 0.0

    call openFile(fileName)
    iostatus = 0
    
    do 
      j = 0
      call readConfiguration(particleArray, radius, height, iostatus)
      if (iostatus < 0) then
        write(*, *) 'End of file ', fileName
        exit
      else if (iostatus > 0) then
        write(*, *) 'Error code ', iostatus
        stop
      else  
        n = size(particleArray)
        write(*, *) 'Particle count ', n
        !! write(*, *) particleArray(1:n)
        do i = 1, size(particleArray)
          if (.not. particleArray(i)%rod) then
            j = j + 1
          end if 
        end do
      end if
      write(*, *) j, ' non-rods in configuration'
    end do

    call closeFile 
  end subroutine testRods


  
  subroutine testNumberOfConfigurations
    use StateReader
    implicit none
    type(particledat), dimension(:), pointer :: particleArray
    real(dp) :: radius, height
    integer :: i
    integer :: iostatus
    character(len = *), parameter :: fileName = 'simdata.out'

    radius = 0.0
    height = 0.0

    call openFile(fileName)
    i = 0
    iostatus = 0

    do 
      call readConfiguration(particleArray, radius, height, iostatus)
      if (iostatus < 0) then
        write(*, *) 'End of file ', fileName
        exit
      else if (iostatus > 0) then
        write(*, *) 'Error code ', iostatus, ' when reading file.'
        stop
      else
        i = i+1
      end if
    end do

    write(*, *) 'There were ', i, ' configurations'
    
    call closeFile
  end subroutine testNumberOfConfigurations



end program TestStateReader

