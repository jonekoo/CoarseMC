
  !Palauttaa tiedostonimen annettujen parametrien perusteella. 
  function filename(R,T,P,GB,Xe,index) result(name)
    implicit none
    intrinsic anint
    character(len=27) :: name
    real(dp),intent(in) :: T,P,R
    integer, intent(in) :: GB,Xe,index
    integer :: Tint,Pint,Rint
    character(len=3) :: Tchar
    character(len=3) :: Rchar
    character(len=7) :: GBchar
    character(len=5) :: Xechar
    character(len=6) :: indexchar
    integer,parameter :: mykind=selected_int_kind(8)
    Tint=anint(100*T,mykind)
    Pint=anint(10*P,mykind)
    Rint=anint(R,mykind)
    write (Tchar,'(I3.3)') Tint
    write (GBchar,'(I7.7)') GB
    write (Xechar,'(I5.5)') Xe
    write (indexchar,'(I6.6)') index
    write (Rchar,'(I3.3)') Rint
    name='R'//Rchar//'T'//Tchar//'.'//trim(adjustl(indexchar))
  end function filename


