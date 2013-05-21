module pov
use nrtype
use particle
use class_poly_box
implicit none

contains

  subroutine povout(simbox, particlearray, povfile)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particlearray
    character(len=*), intent(in) :: povfile

    integer, parameter :: povunit = 18  !! Replace with get_freeunit.
    integer :: opened, N, i
    real(dp) :: x, y, z, ux, uy, uz, camdist, litedist  

    !camdist = maxval((/getx(simbox), gety(simbox), getz(simbox)/))
    litedist = 4._dp * getx(simbox)
    camdist = 3._dp * getx(simbox)
    open(povunit, FILE = povfile, status = 'replace', position = 'append', iostat = opened)
    if(opened /= 0) then
      write(*, *) 'pov:povout:Failed to open ', povfile
      stop
    end if
      write(povunit, *) '#include "colors.inc"'
      write(povunit, *) '#include "transforms.inc"'
      write(povunit, *) 'background {color White}'
      write(povunit, *) 'camera {' 
      write(povunit, *) 'location <',camdist,',',camdist,',',-0,'>'
      write(povunit, *) 'sky <0,1,0>'
      write(povunit, *) 'look_at <0,0,0> }'
      write(povunit, *) 'light_source { <',litedist,',',litedist,',',-litedist,'> color White}'
      !write(povunit, *) 'cylinder{<0,0,',-Lz/2._dp,'>,<0,0,',Lz/2._dp,'>,',R,'open'
      !write(povunit, *) 'texture{pigment{color Yellow}}'
      !write(povunit, *) 'clipped_by{'
      !write(povunit, *) 'plane{x,0}}}'

      N = size(particlearray)
      do i = 1, N
        x = particlearray(i)%x
        y = particlearray(i)%y
        z = particlearray(i)%z
        ux = particlearray(i)%ux
        uy = particlearray(i)%uy
        uz = particlearray(i)%uz
        write(povunit, *) 'sphere {'
        write(povunit, *) '<0,0,0>,0.5'
        if(particlearray(i)%rod) then
          write(povunit, *) 'texture {pigment{color Grey}}'
          write(povunit, *) 'scale <4.4,1,1>'
          write(povunit, *) 'Reorient_Trans(<1,0,0>,<',ux,',',uy,',',-uz,'>)'
        else
          write(povunit, *) 'scale 0.869'
          write(povunit, *) 'texture {pigment{color Red}}'
        end if
        write(povunit, *) 'translate <',x,',',y,',',-z,'>}'
      end do  
      close(povunit)
  end subroutine povout





end module
