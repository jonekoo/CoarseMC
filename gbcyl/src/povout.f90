  !Aliohjelma, joka tuottaa POV-Ray -syöttötiedoston partikkelitaulukosta
  !Tekijä: Jouni Karjalainen
  !R ja Lz määräävät kameran ja valonlähteen koordinaatit. 
  subroutine povout(particlearray, R, Lz)
    implicit none
    type(particledat), dimension(:), pointer :: particlearray
    real(dp), intent(in) :: R, Lz
    integer, parameter :: povunit=18
    character(len=*),parameter :: povfile='povout.pov'
    integer :: opened,N,i
    real(dp) :: x,y,z,ux,uy,uz,camdist,litedist
    camdist=2.0*R
    litedist=camdist
    open(povunit,FILE=povfile,status='replace',position='append',iostat=opened)
    if(opened/=0) then
      write(*,*) 'Tiedoston ',povfile,' avaaminen epäonnistui.'
    else
      write(povunit,*) '#include "colors.inc"'
      write(povunit,*) '#include "transforms.inc"'
      write(povunit,*) 'background {color White}'
      write(povunit,*) 'camera {' 
      write(povunit,*) 'location <',camdist,',',camdist,',',-Lz,'>'
      write(povunit,*) 'sky <0,0,-1>'
      write(povunit,*) 'look_at <0,0,0> }'
      write(povunit,*) 'light_source { <',litedist,',',litedist,',',-Lz,'> color White}'
      write(povunit,*) 'cylinder{<0,0,',-Lz/2,'>,<0,0,',Lz/2,'>,',R,'open'
      write(povunit,*) 'texture{pigment{color Yellow}}'
      write(povunit,*) 'clipped_by{'
      write(povunit,*) 'plane{x,0}}}'

      N=size(particlearray)
      do i=1,N
        x=particlearray(i)%x
        y=particlearray(i)%y
        z=particlearray(i)%z
        ux=particlearray(i)%ux
        uy=particlearray(i)%uy
        uz=particlearray(i)%uz
        write(povunit,*) 'sphere {'
        write(povunit,*) '<0,0,0>,0.5'
        if(particlearray(i)%rod) then
          write(povunit,*) 'texture {pigment{color Grey}}'
          write(povunit,*) 'scale <4.4,1,1>'
          write(povunit,*) 'Reorient_Trans(<1,0,0>,<',ux,',',uy,',',-uz,'>)'
        else
          write(povunit,*) 'scale 0.96'
          write(povunit,*) 'texture {pigment{color Red}}'
        end if
        write(povunit,*) 'translate <',x,',',y,',',-z,'>}'
      end do  

      close(povunit)
    end if
  end subroutine povout



