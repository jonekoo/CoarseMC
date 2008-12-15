module stream_io
  use particle, only : particle_write_state => write_module_state
  use mcstep, only : mcstep_write_state => write_module_state

contains

  subroutine write_state(unit)
    implicit none
    
    integer, intent(in) :: unit
    
    call particle_write_state(unit)
    call mcstep_write_state(unit)
  end subroutine writeState

end module stream_io
