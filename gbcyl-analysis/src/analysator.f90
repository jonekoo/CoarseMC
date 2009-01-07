module analysator

  contains
  
  subroutine analyse(particleArray, np, outputFile)
    use particle, only : particledat
    use order, only : eigens
    implicit none
    
    type(particledat), dimension(np), intent(in) :: particleArray
    integer, intent(in) :: np
    character(len = *), intent(in) :: outputFile
    real(dp), dimesion(3) :: values
    real(dp), dimension(3,3) :: vectors    
    integer, parameter :: writeUnit = 10

    open(writeUnit, FILE=outputFile, STATUS = 'NEW',&
      & ACTION = 'WRITE', POSITION = 'APPEND') 

    call eigens(particleArray, np, values, vectors)
    write(writeUnit, *)  
    
  end subroutine analyse

end module analysator
