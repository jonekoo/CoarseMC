module pov
use num_kind
use m_particledat
use class_poly_box
use utils, only: fmt_char_dp
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
      write(povunit, *) '#version 3.1'
      write(povunit, *) '#include "colors.inc"'
      write(povunit, *) '#include "transforms.inc"'
      write(povunit, *) '#declare GB = sphere {<0, 0, 0>, 0.5 scale <4.4,1,1> texture{pigment{color Grey}}}'
      write(povunit, *) '#declare Xe = sphere {<0, 0, 0>, 0.5 scale 0.869 texture{pigment{color Red}}}'
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
        if(particlearray(i)%rod) then
           write(povunit, fmt='(A)', advance='no') 'object{GB '
           write(povunit, fmt='(A24, 2(' // fmt_char_dp() // ',A1,1X),' // &
                fmt_char_dp() //',A2)', advance='no') &
                'Reorient_Trans(<1,0,0>,<',ux,',',uy,',',-uz,'>)'
        else
           write(povunit, '(A)', advance='no') 'object{Xe '
        end if
        write(povunit, '(A11, 2(' // fmt_char_dp() // ',A1,1X),' // &
             fmt_char_dp() //',A2)') &
             'translate <',x,',',y,',',-z,'>}'
      end do  
      close(povunit)
  end subroutine povout

end module
