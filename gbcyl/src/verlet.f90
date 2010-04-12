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

  public :: initvlist
  public :: updatelist
  public :: freevlist
  public :: pairinteractions
  public :: newlist
  public :: verlet_writeparameters
 
  real(dp), save :: rlist
  integer, save :: nparticles
  integer, save :: nneighbours
  real(dp), dimension(:, :), allocatable, save :: xyzlist
  integer, dimension(:), allocatable, save :: neighbourcounts
  integer, dimension(:, :), allocatable, save ::  neighbours
  integer, save :: nneighboursmax

  interface pairinteractions
    module procedure totpairV, singleparticleV
  end interface

  interface initvlist
    module procedure initvlistparameterizer, initvlistr
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
  end subroutine

  !! Initializes the verlet neighbour list.
  !! 
  !! @p particles the array of particles
  !! @p nparticles the number of particles
  !! @p simbox the simulation cell
  !! @p reader gives the parameters
  !! 
  subroutine initvlistparameterizer(particles, nparticles, simbox, reader)
    type(particledat), dimension(:), intent(in) :: particles
    integer :: nparticles
    type(poly_box), intent(in) :: simbox
    type(parameterizer), intent(in) :: reader
    real(dp) :: rlist
    call getparameter(reader, 'r_list', rlist)
    call initvlist(particles, nparticles, simbox, rlist)
  end subroutine

  !! Initializes the verlet neighbour list. 
  !! 
  !! @p particles the array of particles
  !! @p nparticles the number of particles in @p particles
  !! @p simbox the simulation cell
  !! @p rlist the radius of the neighbourlist
  !! @p rcutoff the cutoff radius of the pair interactions
  !! 
  subroutine initvlistr(particles, nparticlesin, simbox, rlistin)
    intrinsic min
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: nparticlesin
    type(poly_box), intent(in) :: simbox
    real(dp), intent(in) :: rlistin
    integer :: astat
    rlist = rlistin
    nneighboursmax = 500
    nparticles = nparticlesin
    nneighbours = min(nparticles, nneighboursmax) 
    if (nparticles < nneighbours) then
      nneighbours = nparticles
    end if
    allocate(xyzlist(3, nparticles), &
      neighbours(nparticles, nneighbours), & 
      neighbourcounts(nparticles), stat = astat)
    if (astat /= 0) then
      stop 'Could not allocate memory for verlet list.'
    end if
    call newlist(particles, nparticles, simbox)
  end subroutine

  !! Frees the arrays needed by the verlet list.
  !! 
  subroutine freevlist()
    if (allocated(xyzlist)) deallocate(xyzlist)
    if (allocated(neighbours)) deallocate(neighbours)
    if (allocated(neighbourcounts)) deallocate(neighbourcounts)
  end subroutine freevlist

  !! Calculates the total pair interaction energy of all @p particles. Uses
  !! Verlet neighbour list. 
  !!
  !! @p particles the array of all particles
  !! @p nparticles the number of particles
  !! @p simbox the simulation cell
  !! @p Vtot the total pair interaction energy
  !! @p overlap true if some particles overlap with each other
  !! 
  subroutine totpairV(simbox, particles, Vtot, overlap)
    type(particledat), dimension(:), intent(in) :: particles
    type(poly_box), intent(in) :: simbox
    integer :: i, j, jj
    real(dp), intent(out) :: Vtot 
    real(dp) :: pairE
    logical, intent(out) :: overlap
    Vtot = 0._dp
    overlap = .false.
    call updatelist(particles, size(particles), simbox)
    do i = 1, size(particles)
      if(neighbourcounts(i) == 0) cycle
      do jj = 1, neighbourcounts(i)
        j = neighbours(i, jj)
        if(j <= i) cycle
        call pairV(particles(i), particles(j), simbox, pairE, overlap)
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
  subroutine singleparticleV(simbox, particles, particlei, i, singleV, overlap)
    type(particledat), dimension(:), intent(in) :: particles
    type(particledat), intent(in) :: particlei
    integer, intent(in) :: i
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: singleV
    logical, intent(out) :: overlap
    integer :: j
    real(dp) :: pairE
    real(dp) :: distancemoved 
    type(particledat) :: particlevij    
    singleV = 0._dp
    overlap = .false.
    distancemoved = mindistance(simbox, xyzlist(1:3, i), &
    position(particlei))
    if(distancemoved > 0.5_dp * (rlist - cutoff())) then
      call newlist(particles, size(particles), simbox)
    end if
    do j = 1, neighbourcounts(i)
      particlevij = particles(neighbours(i, j))
      call pairV(particlei, particlevij, simbox, pairE, overlap)
      if (overlap) then
        exit
      else
        singleV = singleV + pairE
      end if
    end do
  end subroutine

  !! Makes a new Verlet neighbourlist
  !!
  !! @p particles the array of particles
  !! @p nparticles the number of particles
  !! @p simbox the simulation cell
  !!
  subroutine newlist(particles, nparticles, simbox)
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: nparticles
    type(poly_box), intent(in) :: simbox
    integer :: i, j
    real(dp) :: rij
    do i = 1, nparticles
      neighbourcounts(i) = 0
      xyzlist(1, i) = particles(i)%x
      xyzlist(2, i) = particles(i)%y
      xyzlist(3, i) = particles(i)%z
    end do 
    do i = 1, nparticles - 1
      do j = i + 1, nparticles
        rij = mindistance(simbox, position(particles(i)), &
        position(particles(j)))
        if(rij < rlist ) then
          neighbourcounts(i) = neighbourcounts(i) + 1
          neighbourcounts(j) = neighbourcounts(j) + 1 
          neighbours(i, neighbourcounts(i)) = j
          neighbours(j, neighbourcounts(j)) = i
        end if
      end do
    end do
  end subroutine newlist 

  !! Updates the neighbourlist if needed.
  !!
  !! @p particles the array of particles to be included in the list
  !! @p nparticles the number of particles in @p particles
  !! @p simbox the simulation cell 
  !!
  subroutine updatelist(particles, nparticles, simbox)
    intrinsic max, abs
    type(particledat), dimension(:), intent(in) :: particles     
    integer, intent(in) :: nparticles
    type(poly_box), intent(in) :: simbox
    integer :: i
    do i = 1, nparticles
      if(mindistance(simbox, xyzlist(1:3, i), position(particles(i))) > &
      0.5_dp * (rlist - cutoff())) then
        call newlist(particles, nparticles, simbox)
        exit
      end if
    end do  
  end subroutine updatelist

end module verlet
