module verlet
  !! Module for calculating pair interactions of particles using a Verlet 
  !! neighbour list. 
  use nrtype, only: dp
  use particle, only: particledat, position, getmaxmoves
  use class_poly_box
  use class_parameterizer
  use class_parameter_writer
  use class_pair_potential
  implicit none 
  private  

  public :: verlet_init
  public :: pairinteractions
  public :: verlet_writeparameters
  public :: new_verletlist
  public :: verletlist
  public :: update
  public :: delete
 
  real(dp), save :: rlist = 6.8_dp
  integer, save :: nneighboursmax = 500

  type verletlist
    private
    integer :: nneighbours = 500
    real(dp), dimension(:, :), pointer :: xyzlist => NULL()
    integer, dimension(:), pointer :: neighbourcounts => NULL()
    integer, dimension(:, :), pointer ::  neighbours => NULL()
  end type 

  interface pairinteractions
    module procedure totpairV, singleparticleV
  end interface

  interface verlet_init
    module procedure verlet_initwtparameterizer, verlet_initwtvalues
  end interface

  interface delete
    module procedure verlet_delete
  end interface

  interface update
    module procedure verlet_update
  end interface

  contains

  !! Writes verlet list parameters to parameter file.
  !!
  !! @p writer writes the parameters given by the routine
  !! 
  subroutine verlet_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'Verlet neighbourlist parameters')
    call writeparameter(writer, 'r_list', rlist)
    call writeparameter(writer, 'n_neighbours_max', nneighboursmax)
  end subroutine

  subroutine verlet_initwtparameterizer(reader)
    type(parameterizer), intent(in) :: reader
    call getparameter(reader, 'r_list', rlist)
    call getparameter(reader, 'max_nneighbours', nneighboursmax)
  end subroutine

  subroutine verlet_initwtvalues(rlistin, nneighboursmaxin)
    real(dp), intent(in) :: rlistin
    integer, intent(in) :: nneighboursmaxin
    rlist = rlistin
    nneighboursmax = nneighboursmaxin
  end subroutine

  !! Initializes the verlet neighbour list. 
  !! 
  !! @p particles the array of particles
  !! @p nparticles the number of particles in @p particles
  !! @p simbox the simulation cell
  !! @p rlist the radius of the neighbourlist
  !! @p rcutoff the cutoff radius of the pair interactions
  !! 
  function new_verletlist(simbox, particles) result(vl)
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    type(verletlist) :: vl 
    integer :: nparticles
    integer :: astat
    nparticles = size(particles)
    vl%nneighbours = min(nparticles, nneighboursmax) 
    allocate(vl%xyzlist(3, nparticles), &
      vl%neighbours(nparticles, vl%nneighbours), & 
      vl%neighbourcounts(nparticles), stat = astat)
    if (astat /= 0) then
      stop 'Could not allocate memory for verlet list.'
    end if
    call verlet_update(vl, simbox, particles)
  end function


  !! Makes a new Verlet neighbourlist
  !!
  !! @p particles the array of particles
  !! @p nparticles the number of particles
  !! @p simbox the simulation cell
  !!
  subroutine verlet_update(vl, simbox, particles)
    type(verletlist), intent(inout) :: vl
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    integer :: nparticles
    integer :: i, j
    real(dp) :: rij
    real(dp), dimension(3) :: rijvec
    nparticles = size(particles)
    do i = 1, nparticles
      vl%neighbourcounts(i) = 0
      vl%xyzlist(1, i) = particles(i)%x
      vl%xyzlist(2, i) = particles(i)%y
      vl%xyzlist(3, i) = particles(i)%z
    end do 
    do i = 1, nparticles - 1
      do j = i + 1, nparticles
        !rijvec = minimage(simbox, position(particles(i)), &
        !position(particles(j)))
        rijvec = minimage(simbox, position(particles(j)) - &
        position(particles(i)))
        rij = sqrt(dot_product(rijvec, rijvec))
        if(rij < rlist ) then
          vl%neighbourcounts(i) = vl%neighbourcounts(i) + 1
          vl%neighbourcounts(j) = vl%neighbourcounts(j) + 1 
          vl%neighbours(i, vl%neighbourcounts(i)) = j
          vl%neighbours(j, vl%neighbourcounts(j)) = i
        end if
      end do
    end do
  end subroutine 

  subroutine verlet_delete(vl)
    type(verletlist), intent(inout) :: vl
    if (associated(vl%xyzlist)) deallocate(vl%xyzlist)
    if (associated(vl%neighbours)) deallocate(vl%neighbours)
    if (associated(vl%neighbourcounts)) deallocate(vl%neighbourcounts)
  end subroutine

  !! Calculates the total pair interaction energy of all @p particles. Uses
  !! Verlet neighbour list. 
  !!
  !! @p particles the array of all particles
  !! @p nparticles the number of particles
  !! @p simbox the simulation cell
  !! @p Vtot the total pair interaction energy
  !! @p overlap true if some particles overlap with each other
  !! 
  subroutine totpairV(vl, simbox, particles, Vtot, overlap)
    type(verletlist), intent(in) :: vl
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_box), intent(in) :: simbox
    integer :: i, j, jj
    real(dp), intent(out) :: Vtot 
    real(dp) :: pairE
    logical, intent(out) :: overlap
    Vtot = 0._dp
    overlap = .false.
    do i = 1, size(particles)
      if(vl%neighbourcounts(i) == 0) cycle
      do jj = 1, vl%neighbourcounts(i)
        j = vl%neighbours(i, jj)
        if(j <= i) cycle
        call pairv(particles(i), particles(j), simbox, pairE, overlap)
        if (overlap) return
        Vtot = Vtot + pairE      
      end do
    end do
  end subroutine

  !! Calculates the interaction energy of @p particlei with all the other 
  !! particles.
  !!
  !! @p particles the array of particles
  !! @p nparticles the number of particles in @p particles
  !! @p simbox the simulation cell
  !! @p particlei the particle to which the interactions are calculated
  !! @p i the index of @p particlei in @p particles
  !! @p singleV the interaction energy of @p particlei with other particles.
  !! @p overlap true if particlei and some other particle overlap
  !!
  subroutine singleparticleV(vl, simbox, particles, i, singleV, overlap)
    type(verletlist), intent(in) :: vl
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: i
    real(dp), intent(out) :: singleV
    logical, intent(out) :: overlap
    integer :: j
    real(dp) :: pairE
    type(particledat) :: particlevij    
    singleV = 0._dp
    overlap = .false.
    do j = 1, vl%neighbourcounts(i)
      particlevij = particles(vl%neighbours(i, j))
      call pairV(particles(i), particlevij, simbox, pairE, overlap)
      if (overlap) then
        exit
      else
        singleV = singleV + pairE
      end if
    end do
  end subroutine

end module verlet
