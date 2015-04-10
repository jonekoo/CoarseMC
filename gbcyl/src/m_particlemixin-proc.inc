
  subroutine MIXINPARAM_energy(this, res, err)
    class(particlemixin), intent(inout) :: this
    real(REAL64), intent(out) :: res
    integer, intent(out) :: err
    real(REAL64) :: temp
    integer :: i
    res = 0
    if (associated(this%interactions)) then
       do i = 1, size(this%interactions)
          call this%interactions(i)%value(this, temp, err)
          if (err /= 0) return
          res = res + temp
       end do
    else
       write(error_unit, *) 'Error: No interactions defined for ' // this%name
    end if
  end subroutine MIXINPARAM_energy

  subroutine particlemixin_setname(name)
    type(particlemixin) :: particlemixin_name
    character(len=*), intent(in) :: name
    particlemixin_name%name = name
  end subroutine
