module energy
  use nrtype, only: dp
  use particle
  use particlewall
  use class_simplelist
  use class_poly_box
  use class_parameter_writer
  use class_parameterizer
  use class_pair_potential
  implicit none
  private
 
  public :: potentialenergy
  public :: totwallprtclV
  public :: energy_init
  public :: energy_writeparameters
  public :: get_cutoff

  logical, save :: iswall
  real(dp), save :: rcutoff = 5.5_dp
  logical :: is_initialized = .false.

  interface potentialenergy
    module procedure totalenergy, singleenergy, simple_singleenergy
  end interface

  contains

  real(dp) function get_cutoff()
    if (.not. is_initialized) stop 'Trying to access all_pairs:get_cutoff before the module is initialized.'
    get_cutoff = rcutoff
  end function

  subroutine energy_init(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'is_wall_on', iswall)
    if (iswall) call initptwall(reader)
    call getparameter(reader, 'r_cutoff', rcutoff)
    call pp_init(reader)
    is_initialized = .true.
  end subroutine
 
  subroutine energy_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writeparameter(writer, 'is_wall_on', iswall)
    if (iswall) call particlewall_writeparameters(writer)
    call writeparameter(writer, 'r_cutoff', rcutoff)
    call pp_writeparameters(writer)
  end subroutine

  !> Calculates the total potential energy of @param particles with respect to 
  !! each other and any external fields. 
  !!
  !! @param simbox the simulation box where the particles are.
  !! @param particles the array of particles to which the energy is calculated.
  !! @param nbrlist the neighbourlist used to calculate interactions of 
  !! particles with each other. 
  !! @param Etot the total energy when the routine returns.
  !! @param overlap tells if there's an overlap meaning an infinite 
  !! interaction energy.
  !! 
  subroutine totalenergy(simbox, particles, Etot, overlap, n_pairs)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    real(dp), intent(out) :: Etot
    logical, intent(out) :: overlap
    integer, intent(out), optional :: n_pairs
    real(dp) :: Vpairtot
    real(dp) :: Vwalltot
    Etot = 0._dp
    Vwalltot = 0._dp
    Vpairtot = 0._dp
    overlap = .false.
    call allpairinteractions(simbox, particles, Vpairtot, overlap, n_pairs)
    if (iswall .and. (.not. overlap)) then
      call totwallprtclV(simbox, particles, Vwalltot, overlap)
    end if  
    if (.not. overlap) then
      Etot = Vpairtot + Vwalltot
    end if
  end subroutine totalenergy
          
  !> Accumulates the particle-wall interaction for @param particles.
  !! 
  !! @param simbox the simulation box where the particles are.
  !! @param particles the array of particles.
  !! @param Eptwlltot the total particle-wall interaction at return.
  !! @param overlap is true at return only if the interaction is infinitely 
  !! large.
  !! 
  subroutine totwallprtclV(simbox, particles, Eptwlltot, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    real(dp), intent(out) :: Eptwlltot
    logical, intent(out) :: overlap 
    integer :: i
    real(dp) :: oneprtclV 
    Eptwlltot = 0._dp
    overlap = .false.
    !! Could this be replaced with forall or where construct?
    !! If volume scaling trials are done only in the direction of the axis of
    !! the cylinder, then yes since no overlap can occur. It is also possible
    !! otherwise but may not give any performance benefit. 
    do i = 1, size(particles)
      call particlewall_potential(particles(i), simbox, oneprtclV, overlap)
      if (overlap) return
      Eptwlltot = Eptwlltot + oneprtclV
    end do
  end subroutine

  !> Returns the potential energy of @param particles(@param i) in the system.
  !! 
  !! @param simbox the simulation box where the particles are.
  !! @param particles the particles in the system.
  !! @param nbrlist the neighbourlist to be used to calculate short-ranged
  !! pair interactions.
  !! @param i the index of the particle in @param particles.
  !! @param Vitot the potential energy at return.
  !! @param overlap is .true. if an infinitely large interaction occurs. 
  !!
  pure subroutine singleenergy(simbox, particles, nbrlist, i, Vitot, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    type(simplelist), intent(in) :: nbrlist
    integer, intent(in) :: i
    real(dp), intent(out) :: Vitot
    logical, intent(out) :: overlap
    real(dp) :: Vipair 
    real(dp) :: Viwall 
    Viwall = 0._dp 
    Vipair = 0._dp
    Vitot = 0._dp
    overlap = .false.
    if (iswall) then
      call particlewall_potential(particles(i), simbox, Viwall, overlap)
    end if
    if (.not. overlap) then
      call pairinteractions(simbox, particles, i, Vipair, overlap)
    end if
    if (.not. overlap) then 
      Vitot = Vipair + Viwall
    end if
  end subroutine

  !> Returns the potential energy of @param particles(@param i) in the system.
  !! 
  !! @param simbox the simulation box where the particles are.
  !! @param particles the particles in the system.
  !! @param nbrlist the neighbourlist to be used to calculate short-ranged
  !! pair interactions.
  !! @param i the index of the particle in @param particles.
  !! @param Vitot the potential energy at return.
  !! @param overlap is .true. if an infinitely large interaction occurs. 
  !!
  pure subroutine simple_singleenergy(simbox, particles, i, Vitot, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: i
    real(dp), intent(out) :: Vitot
    logical, intent(out) :: overlap
    real(dp) :: Vipair 
    real(dp) :: Viwall 
    Viwall = 0._dp 
    Vipair = 0._dp
    Vitot = 0._dp
    overlap = .false.
    if (iswall) then
      call particlewall_potential(particles(i), simbox, Viwall, overlap)
    end if
    if (.not. overlap) then
      call pairinteractions(simbox, particles, i, Vipair, overlap)
    end if
    if (.not. overlap) then 
      Vitot = Vipair + Viwall
    end if
  end subroutine

  pure subroutine allpairinteractions(simbox, particles, energy, overlap, n_pairs)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    integer, intent(out), optional :: n_pairs
    integer :: nparticles
    integer :: i, j
    real(dp) :: epair
    real(dp) :: rij(3)
    if (present(n_pairs)) n_pairs = 0 
    energy = 0._dp
    overlap = .false.
    nparticles = size(particles)
    do i = 1, nparticles - 1
      do j = i + 1, nparticles
        rij = minimage(simbox, position(particles(j))-position(particles(i)))
        if (dot_product(rij, rij) < rcutoff**2) then
          call pairV(particles(i), particles(j), rij, epair, overlap)
          if (present(n_pairs)) n_pairs = n_pairs + 1
          if(overlap) then
            return
          else
            energy = energy + epair
          end if
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
  pure subroutine pairinteractions(simbox, particles, i, energy, overlap, n_pairs)
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_box), intent(in) :: simbox
    integer, intent(in) :: i
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    integer, intent(out), optional :: n_pairs
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
    if (present(n_pairs)) n_pairs = 0

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
        call pairV(particles(i), particles(j), rijs(:,j), epair, &
          overlap)
        if (present(n_pairs)) n_pairs = n_pairs + 1
        !overlap = overlap .or. overlap_ij
        if(overlap) then
          !write(*, *) 'Overlap in all_pairs:singleparticlepairs'
          return
        else
          energy = energy + epair
        end if
      end if
    end do
  end subroutine

end module
