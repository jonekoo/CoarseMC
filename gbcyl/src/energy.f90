module energy
  use nrtype, only: dp
  use particle, only: particledat
  use particlewall
  use class_poly_nbrlist
  use class_poly_box
  use class_parameter_writer
  use class_parameterizer
  implicit none
  private
 
  public :: potentialenergy
  public :: totwallprtclV
  public :: energy_init
  public :: energy_writeparameters

  logical, save :: iswall

  interface potentialenergy
    module procedure totalenergy, singleenergy
  end interface

  contains

  subroutine energy_init(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'is_wall_on', iswall)
    if (iswall) call initptwall(reader)
  end subroutine
 
  subroutine energy_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writeparameter(writer, 'is_wall_on', iswall)
    if (iswall) call particlewall_writeparameters(writer)
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
  subroutine totalenergy(simbox, particles, nbrlist, Etot, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_nbrlist), intent(in) :: nbrlist
    real(dp), intent(out) :: Etot
    logical, intent(out) :: overlap
    real(dp) :: Vpairtot
    real(dp) :: Vwalltot
    Etot = 0._dp
    Vwalltot = 0._dp
    overlap = .false.
    call pairinteractions(nbrlist, simbox, particles, Vpairtot, overlap)
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
  subroutine singleenergy(simbox, particles, nbrlist, i, Vitot, overlap)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_nbrlist), intent(in) :: nbrlist
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
      call pairinteractions(nbrlist, simbox, particles, i, Vipair, overlap)
    end if
    if (.not. overlap) then 
      Vitot = Vipair + Viwall
    end if
  end subroutine

end module
