  !Aliohjelma, joka kirjoittaa partikkelien tiedot tiedostoon 
  !Tekijä: Jouni Karjalainen 
  !Muokattu Juho Lintuvuoren alkuperäisestä koodista. 
  subroutine writeatomsout(N,radius,height,particle)
    implicit none
    type(particledat), dimension(:), pointer :: particle
    real(dp),intent(in) :: radius,height
    integer :: N
    integer :: ios,i
    integer,parameter :: wunit=20
    type(particledat), pointer :: particlei
    character(len=40), parameter :: wfile='atoms.out'  

    open(wunit,file=wfile,status='replace',form='formatted',iostat=ios)
    if(ios/=0)then
       write(*,*)'Virhe tiedoston ',wfile,'avaamisessa. Ohjelma keskeytyy..'
       stop;
    end if
    N=size(particle)      
    write(wunit,*)N,radius,height
    do i=1,N
      particlei=>particle(i)
      write(wunit,*) particlei
    end do
    close(wunit);
  end subroutine writeatomsout


