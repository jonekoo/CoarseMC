module all_pairs
  use particle
  use nrtype
  use class_poly_box
  use class_pair_potential
  use class_parameterizer
  use class_parameter_writer
  implicit none
  private

  public :: pairinteractions
  public :: stublist
  public :: ap_init
  public :: ap_writeparameters
  public :: cutoff

  type stublist
    private
    character(len = 30) :: name = 'stublist'
  end type

  interface pairinteractions
    module procedure singleparticlepairs, ap_allpairinteractions, ap_singleparticlepairs!, allpairinteractions
  end interface

  real(dp), save :: rcutoff = 5.5_dp

  contains

  elemental function cutoff() result(r)
    real(dp) :: r
    r = rcutoff
  end function

  subroutine ap_init(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'r_cutoff', rcutoff)
  end subroutine

  subroutine ap_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writeparameter(writer, 'r_cutoff', rcutoff)
  end subroutine

  subroutine ap_allpairinteractions(sl, simbox, particles, energy, overlap)
    type(stublist), intent(in) :: sl
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    call allpairinteractions(simbox, particles, energy, overlap)
  end subroutine

  pure subroutine ap_singleparticlepairs(sl, simbox, particles, i, energy, overlap)
    type(stublist), intent(in) :: sl
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: i
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    call pairinteractions(simbox, particles, i, energy, overlap)
  end subroutine

  pure subroutine allpairinteractions(simbox, particles, energy, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    integer :: nparticles
    integer :: i, j
    real(dp) :: epair
    real(dp) :: rij(3)
    energy = 0._dp
    overlap = .false.
    nparticles = size(particles)
    do i = 1, nparticles - 1
      do j = i + 1, nparticles
        rij = minimage(simbox, position(particles(j))-position(particles(i)))
        call pairV(particles(i), particles(j), rij, simbox, epair, overlap)
        if(overlap) then
          return
        else
          energy = energy + epair
        end if
      end do
    end do
  end subroutine

  !> Calls the pair potential calculation for particles(i) and all particles in
  !! @p particles for which the center-to-center distance is less than cutoff
  !! radius.
  !! 
  !! @param simbox the simulation box for the minimum image calculation.
  !! @param particles the array of particles.
  !! @param i the particle to which the interactions are calculated.
  !! @param energy sum of all pair interactions for particles(i) and other
  !! particles.
  !! @param overlap is true if some two particles overlap. Correct energy is not
  !! guaranteed when overlap = .true. 
  !! 
  pure subroutine singleparticlepairs(simbox, particles, i, energy, overlap)
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_box), intent(in) :: simbox
    integer, intent(in) :: i
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    integer :: nparticles
    integer :: j
    real(dp) :: epair
    logical :: cutoff_mask(size(particles))
    real(dp) :: rijs(3, size(particles))
    logical :: overlap_ij
    energy = 0._dp
    overlap = .false.
    overlap_ij = .false.
    nparticles = size(particles)

    !! Select the calculated interactions by cutoff already here
    do j = 1, nparticles
      !! :NOTE: The chosen way is better than
      !! :NOTE: rij = minimage(simbox, position(particlei), position(particlej))
      !! It is notably faster to use direct references to particle coordinates.
      rijs(:, j) = minimage(simbox, (/particles(j)%x-particles(i)%x,&
        particles(j)%y-particles(i)%y, particles(j)%z-particles(i)%z/))
    end do
    do j = 1, nparticles
      cutoff_mask(j) = dot_product(rijs(:,j), rijs(:,j)) < rcutoff**2
    end do
    !! Remove the particle i from the cutoff mask 
    cutoff_mask(i) = .false. 

    !! :TODO: test overlap conditions also here?

    !! :TODO: sort particles by type here as well?  
    do j = 1, nparticles
      if (cutoff_mask(j)) then
      !if (dot_product(rijs(:, j), rijs(:, j)) < rcutoff**2 .and. j /= i) then
        call pairV(particles(i), particles(j), rijs(:,j), simbox, epair, &
          overlap)
        overlap = overlap .or. overlap_ij
        !if(overlap) then
        !  return
        !else
          energy = energy + epair
        !end if
      end if
    end do
  end subroutine

end module
