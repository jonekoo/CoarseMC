module all_pairs
  use particle
  use nrtype
  use class_poly_box
  use class_pair_potential
  implicit none
  private

  public :: pairinteractions
  public :: stublist
  public :: new_stublist

  type stublist
    private
    character(len = 30) :: name = 'stublist'
  end type

  interface pairinteractions
    module procedure allpairinteractions, singleparticlepairs, sl_allpairinteractions, sl_singleparticlepairs
  end interface

  contains

  function new_stublist()
    type(stublist) :: new_stublist
  end function

  subroutine sl_allpairinteractions(sl, simbox, particles, energy, overlap)
    type(stublist), intent(in) :: sl
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    call pairinteractions(simbox, particles, energy, overlap)
  end subroutine

  subroutine sl_singleparticlepairs(sl, simbox, particles, i, energy, overlap)
    type(stublist), intent(in) :: sl
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: i
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    call pairinteractions(simbox, particles, i, energy, overlap)
  end subroutine

  subroutine allpairinteractions(simbox, particles, energy, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    integer :: nparticles
    integer :: i, j
    real(dp) :: epair
    energy = 0._dp
    overlap = .false.
    nparticles = size(particles)
    do i = 1, nparticles - 1
      do j = i + 1, nparticles
        call pairV(particles(i), particles(j), simbox, epair, overlap)
        if(overlap) then
          return
        else
          if(epair < -1000._dp) write(*, *) 'epair = ', epair
          energy = energy + epair
        end if
      end do
    end do
  end subroutine

  subroutine singleparticlepairs(simbox, particles, i, energy, overlap) 
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_box), intent(in) :: simbox
    integer, intent(in) :: i
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    integer :: nparticles
    integer :: j
    real(dp) :: epair
    energy = 0._dp
    overlap = .false.
    nparticles = size(particles)
    do j = 1, nparticles
      if (j /= i) then
        call pairV(particles(i), particles(j), simbox, epair, overlap)
        if(overlap) then
          return
        else
          energy = energy + epair
        end if
      end if
    end do
  end subroutine

end module
