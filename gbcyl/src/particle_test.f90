module particle_test
use ftnunit
use particle
use num_kind
use mtmod
implicit none

public :: test_all

contains

  subroutine test_all
    call test(test_writeread, "Write a GB particle to file and read it back.")
  end subroutine

  subroutine test_writeread()
    type(particledat) :: particlewritten
    type(particledat) :: particleread
    integer :: unit = 10
    real(dp) :: margin = 1.e-9
    character(len=200) :: line
    integer :: ios
    !! Set random values to the attributes of the particle.
    particlewritten%x = 100._dp * (2._dp * grnd() - 1._dp)
    particlewritten%y = 100._dp * (2._dp * grnd() - 1._dp)
    particlewritten%z = 100._dp * (2._dp * grnd() - 1._dp)
    particlewritten%ux = 2._dp * grnd() - 1._dp
    particlewritten%uy = 2._dp * grnd() - 1._dp
    particlewritten%uz = 2._dp * grnd() - 1._dp
    !! Write to a scratch file.
    open(unit, status='SCRATCH')
    call writeparticle(unit, particlewritten)
    !! Read back from the scratch file.
    rewind unit
    call readparticle(unit, particleread, ios)
    call assert_equal(0, ios, "There was an error reading the particle.")
    !! Compare components with the original
    call assert_comparable(particlewritten%x, particleread%x, &
    margin, "x-coordinate values are not comparable")
    call assert_comparable(particlewritten%y, particleread%y, &
    margin, "x-coordinate values are not comparable")
    call assert_comparable(particlewritten%z, particleread%z, &
    margin, "x-coordinate values are not comparable")
    call assert_comparable(particlewritten%ux, particleread%ux, &
    margin, "x-coordinate values are not comparable")
    call assert_comparable(particlewritten%uy, particleread%uy, &
    margin, "x-coordinate values are not comparable")
    call assert_comparable(particlewritten%uz, particleread%uz, &
    margin, "x-coordinate values are not comparable")
    close(unit)
  end subroutine test_writeread

end module
