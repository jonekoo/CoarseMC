module nrutil
  use nrtype
  implicit none
  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8 !Each NPAR2 must <= 
  INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2  !the corresponding NPAR.
  INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
  INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
  INTEGER(I4B), PARAMETER :: NPAR_POLY=8
  INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
  
  !Next, generic interfaces for routines with overloaded versions.
  !Naming conventions for appended codes in the names of overloaded 
  !routines are as follows: r=real, d=double precision, i=integer, 
  !c=complex, z=double-precision complex, h=character, l=logical. Any 
  !of r,d,i,c,z,h,l may be followed by v=vector or m=matrix (v,m 
  !suffixes are used only when needed to resolve ambiguities)

  INTERFACE swap
    MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
        swap_cv,swap_cm,&!swap_z,swap_zv,swap_zm, &
        masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE
  
  
  interface reallocate
    module procedure reallocate_rv,reallocate_rm,&
      reallocate_iv,reallocate_im,reallocate_hv
  end interface

  interface assert
    MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  end interface
    
  interface assert_eq
    module procedure assert_eq2,assert_eq3,assert_eq4,assert_eqn
  end interface  

  interface arth
    module procedure arth_r, arth_i !, arth_d
  end interface

  interface cumsum
    module procedure cumsum_r, cumsum_i
  end interface

  INTERFACE outerprod
    MODULE PROCEDURE outerprod_r !,outerprod_d
  END INTERFACE

  INTERFACE outerdiff
    MODULE PROCEDURE outerdiff_r,outerdiff_i !,outerdiff_d
  END INTERFACE

contains

  SUBROUTINE swap_i(a,b)
  !Swap the contents of a and b.
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
   
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i


  SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r



  SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum

    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv



  SUBROUTINE swap_c(a,b)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
  
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_c




  SUBROUTINE swap_cv(a,b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum

    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cv



  SUBROUTINE swap_cm(a,b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
    
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cm



  SUBROUTINE swap_z(a,b)
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    COMPLEX(DPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z




  SUBROUTINE swap_zv(a,b)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum

    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv



  SUBROUTINE swap_zm(a,b)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum

    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm



  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(SP) :: swp

    if (mask) then
      swp=a
      a=b
      b=swp
    end if
  END SUBROUTINE masked_swap_rs



  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a)) :: swp

    where (mask)
      swp=a
      a=b
      b=swp
    end where
  END SUBROUTINE masked_swap_rv




  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp

    where (mask)
      swp=a
      a=b
      b=swp
    end where
  END SUBROUTINE masked_swap_rm



  function reallocate_rv(p,n)
    !Reallocate a pointer to a new size, preserving its previous contents. 
    real(sp), dimension(:), pointer :: p, reallocate_rv
    integer(I4B), intent(in) :: n
    integer(I4B) :: nold,ierr
    
    allocate(reallocate_rv(n),stat=ierr)
    if (ierr /= 0) call &
        nrerror('reallocate_rv: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p)
    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  end function reallocate_rv


  function reallocate_iv(p,n)
    integer(I4B), dimension(:), pointer :: p,reallocate_iv
    integer(I4B), intent(in) :: n
    integer(I4B) :: nold, ierr
    
    allocate(reallocate_iv(n), stat=ierr)
    if (ierr /= 0) call &
        nrerror('reallocate_iv: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p)
    reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  end function reallocate_iv


  function reallocate_hv(p,n)
    character(1), dimension(:), pointer :: p, reallocate_hv
    integer(I4B), intent(in) :: n
    integer(I4B) :: nold,ierr
    
    allocate(reallocate_hv(n), stat=ierr)
    if (ierr /= 0) call &
        nrerror('reallocate_hv: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p)
    reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  end function reallocate_hv


  function reallocate_rm(p,n,m)
    real(sp), dimension(:,:), pointer :: p, reallocate_rm
    integer(I4B), intent(in) :: n,m
    integer(I4B) :: nold,mold,ierr
    
    allocate(reallocate_rm(n,m), stat=ierr)
    if (ierr /= 0) call &
        nrerror('reallocate_rm: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p,1)
    mold=size(p,2) 
    reallocate_rm(1:min(nold,n),1:min(mold,m))=&
        p(1:min(nold,n),1:min(mold,m)) 
    deallocate(p)
  end function reallocate_rm


  function reallocate_im(p,n,m)
    integer(I4B), dimension(:,:), pointer :: p, reallocate_im
    integer(I4B), intent(in) :: n,m
    integer(I4B) :: nold,mold,ierr
    
    allocate(reallocate_im(n,m), stat=ierr)
    if (ierr /= 0) call &
        nrerror('reallocate_im: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p,1)
    mold=size(p,2) 
    reallocate_im(1:min(nold,n),1:min(mold,m))=&
        p(1:min(nold,n),1:min(mold,m)) 
    deallocate(p)
  end function reallocate_im


  FUNCTION iminloc(arr)
  !Index of minloc on an array.
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(1) :: imin
    INTEGER(I4B) :: iminloc

    imin=minloc(arr(:))
    iminloc=imin(1)
  END FUNCTION iminloc

  !! Routines for argument checking and error handling:

  subroutine assert1(n1,string)
    !! Report and die if any logical is false (used for arg range checking).
    character(len=*), intent(in) :: string
    logical, intent(in) :: n1
    if (.not. n1) then
      write (*,*) 'nrerror: an assetion failed with this tag:', string
      STOP 'program terminated by assert1'
    end if
  end subroutine assert1

  subroutine assert2(n1,n2,string)
    character(len=*), intent(in) :: string
    logical, intent(in) :: n1,n2
    if (.not. (n1 .and. n2)) then
      write (*,*) 'nrerror: an assertion failed with this tag:', string
      STOP 'program terminated by assert2'
    end if
  end subroutine assert2

  subroutine assert3(n1,n2,n3,string)
    character(len=*), intent(in) :: string
    logical, intent(in) :: n1,n2,n3
    if(.not. (n1 .and. n2 .and. n3)) then
      write(*,*) 'nrerror: an assertion failed with this tag:',string
      STOP 'program terminated by assert3'
    end if
  end subroutine assert3

  subroutine assert4(n1,n2,n3,n4,string)
    character(len=*), intent(in) :: string
    logical, intent(in) :: n1,n2,n3,n4
    if(.not. (n1 .and. n2 .and. n3 .and. n4)) then
      write(*,*) 'nrerror: an assertion failed with this tag:',string
      STOP 'program terminated by assert3'
    end if
  end subroutine assert4

  subroutine assert_v(n,string)
    character(len=*), intent(in) :: string
    logical, dimension(:), intent(in) :: n
    if(.not. all(n)) then
      write(*,*) 'nrerror: an assertion failed with this tag:',string
      STOP 'program terminated by assert_v'
    end if
  end subroutine assert_v



  function assert_eq2(n1,n2,string)
    !Report and die if integers not all equal (used for size checking)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2
    integer :: assert_eq2
    
    if (n1==n2) then
      assert_eq2=n1
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag:',string
      stop 'program terminated by assert_eq2'
    end if
  end function assert_eq2


  function assert_eq3(n1,n2,n3,string)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2,n3
    integer :: assert_eq3
    
    if (n1==n2 .and. n2==n3) then
      assert_eq3=n1
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag:',string
      stop 'program terminated by assert_eq3'
    end if
  end function assert_eq3


  function assert_eq4(n1,n2,n3,n4,string)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2,n3,n4
    integer :: assert_eq4
    
    if (n1==n2 .and. n2==n3 .and. n3==n4) then
      assert_eq4=n1
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag:',string
      stop 'program terminated by assert_eq4'
    end if
  end function assert_eq4


  
  function assert_eqn(nn,string)
    character(len=*), intent(in) :: string
    integer, dimension(:), intent(in) :: nn
    integer :: assert_eqn
    
    if (all(nn(2:)==nn(1))) then
      assert_eqn=nn(1)
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag:',string
      stop 'program terminated by assert_eqn'
    end if
  end function assert_eqn



  subroutine nrerror(string)
    !Reports a messages, then dies. 
    character(len=*), intent(in) :: string
    write(*,*) 'nrerror: ',string
    stop 'program terminated by nrerror'
  end subroutine nrerror



  FUNCTION arth_r(first,increment,n)
!Array function returning an arithmetic progression.
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4B) :: k,k2
    REAL(SP) :: temp

    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
      do k=2,n
        arth_r(k)=arth_r(k-1)+increment
      end do
    else
      do k=2,NPAR2_ARTH
        arth_r(k)=arth_r(k-1)+increment
      end do
      temp=increment*NPAR2_ARTH
      k=NPAR2_ARTH
      do
        if (k >= n) exit
        k2=k+k
        arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
        temp=temp+temp
        k=k2
      end do
    end if
  END FUNCTION arth_r

  

  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
      do k=2,n
        arth_d(k)=arth_d(k-1)+increment
      end do
    else
      do k=2,NPAR2_ARTH
        arth_d(k)=arth_d(k-1)+increment
      end do
      temp=increment*NPAR2_ARTH
      k=NPAR2_ARTH
      do
        if (k >= n) exit
        k2=k+k
        arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
        temp=temp+temp
        k=k2
      end do
    end if
END FUNCTION arth_d



FUNCTION arth_i(first,increment,n)
  INTEGER(I4B), INTENT(IN) :: first,increment,n
  INTEGER(I4B), DIMENSION(n) :: arth_i
  INTEGER(I4B) :: k,k2,temp

  !write(*,*) 'Suoritetaan funktiota arth_i'
  if (n > 0) arth_i(1)=first
  if (n <= NPAR_ARTH) then
    do k=2,n
      arth_i(k)=arth_i(k-1)+increment
    end do
  else
    do k=2,NPAR2_ARTH
      arth_i(k)=arth_i(k-1)+increment
    end do
    temp=increment*NPAR2_ARTH
    k=NPAR2_ARTH
    do
      if (k >= n) exit
      k2=k+k
      arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
      temp=temp+temp
      k=k2
    end do
  end if
  !write(*,*) 'arth_i:poistuttiin funktiosta'
END FUNCTION arth_i



RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
!Cumulative sum on an array, with optional additive seed.
  REAL(SP), DIMENSION(:), INTENT(IN) :: arr
  REAL(SP), OPTIONAL, INTENT(IN) :: seed
  REAL(SP), DIMENSION(size(arr)) :: ans
  INTEGER(I4B) :: n,j
  REAL(SP) :: sd

  n=size(arr)
  if (n == 0_i4b) RETURN
  sd=0.0_sp
  if (present(seed)) sd=seed
  ans(1)=arr(1)+sd
  if (n < NPAR_CUMSUM) then
    do j=2,n
      ans(j)=ans(j-1)+arr(j)
    end do
  else
    ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
    ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
  end if
END FUNCTION cumsum_r


RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
  INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
  INTEGER(I4B), DIMENSION(size(arr)) :: ans
  INTEGER(I4B) :: n,j,sd

  n=size(arr)
  if (n == 0_i4b) RETURN
  sd=0_i4b
  if (present(seed)) sd=seed
  ans(1)=arr(1)+sd
  if (n < NPAR_CUMSUM) then
    do j=2,n
      ans(j)=ans(j-1)+arr(j)
    end do
  else
    ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
    ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
  end if
END FUNCTION cumsum_i


FUNCTION outerprod_r(a,b)
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
  REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r

  outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
  spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprod_r


FUNCTION outerprod_d(a,b)
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
  REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
  
  outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
  spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprod_d



FUNCTION outerdiff_r(a,b)
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
  REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r

  outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
  spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiff_r



FUNCTION outerdiff_d(a,b)
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
  REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d

  outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
  spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiff_d



FUNCTION outerdiff_i(a,b)
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
  INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i

  outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
  spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiff_i



FUNCTION upper_triangle(j,k,extra)
!Return an upper triangular logical mask.
  INTEGER(I4B), INTENT(IN) :: j,k
  INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
  LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
  INTEGER(I4B) :: n

  !write(*,*) 'Suoritetaan funktiota upper_triangle'
  n=0
  if (present(extra)) n=extra
  upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
  !write(*,*) 'Poistutaan funktiosta upper_triangle'
END FUNCTION upper_triangle


end module nrutil
