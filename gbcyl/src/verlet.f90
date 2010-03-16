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
  public :: pair_interactions
  public :: newlist
  public :: verlet_write_parameters
 
  real(dp), save :: r_list_
  integer, save :: n_particles_
  integer, save :: n_neighbours_
  real(dp), dimension(:, :), allocatable, save :: xyz_list_
  integer, dimension(:), allocatable, save :: neighbour_counts_
  integer, dimension(:, :), allocatable, save ::  neighbours_
  integer, save :: n_neighbours_max_

  interface pair_interactions
    module procedure totpairV, singleparticleV
  end interface

  interface initvlist
    module procedure initvlist_parameterizer, initvlist_r
  end interface

  contains

  !! Writes verlet list parameters to parameter file.
  !!
  !! @p writer writes the parameters given by the routine
  !! 
  subroutine verlet_write_parameters(writer)
    type(parameter_writer), intent(in) :: writer
    call write_comment(writer, 'Verlet neighbourlist parameters')
    call write_parameter(writer, 'r_list', r_list_)
  end subroutine

  !! Initializes the verlet neighbour list.
  !! 
  !! @p particles the array of particles
  !! @p n_particles the number of particles
  !! @p simbox the simulation cell
  !! @p reader gives the parameters
  !! 
  subroutine initvlist_parameterizer(particles, n_particles, simbox, reader)
    type(particledat), dimension(:), intent(in) :: particles
    integer :: n_particles
    type(poly_box), intent(in) :: simbox
    type(parameterizer), intent(in) :: reader
    real(dp) :: r_list
    real(dp) :: r_cutoff
    call get_parameter(reader, 'r_list', r_list)
    call initvlist(particles, n_particles, simbox, r_list)
  end subroutine

  !! Initializes the verlet neighbour list. 
  !! 
  !! @p particles the array of particles
  !! @p n_particles the number of particles in @p particles
  !! @p simbox the simulation cell
  !! @p r_list the radius of the neighbourlist
  !! @p r_cutoff the cutoff radius of the pair interactions
  !! 
  subroutine initvlist_r(particles, n_particles, simbox, r_list)
    intrinsic min
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    type(poly_box), intent(in) :: simbox
    real(dp), intent(in) :: r_list
    integer :: astat
    r_list_ = r_list
    n_neighbours_max_ = 500
    n_particles_ = n_particles
    n_neighbours_ = min(n_particles_, n_neighbours_max_) 
    if (n_particles_ < n_neighbours_) then
      n_neighbours_ = n_particles_
    end if
    allocate(xyz_list_(3, n_particles_), &
      neighbours_(n_particles_, n_neighbours_), & 
      neighbour_counts_(n_particles_), stat = astat)
    if (astat /= 0) then
      stop 'Could not allocate memory for verlet list.'
    end if
    call newlist(particles, n_particles, simbox)
  end subroutine

  !! Frees the arrays needed by the verlet list.
  !! 
  subroutine freevlist()
    if (allocated(xyz_list_)) deallocate(xyz_list_)
    if (allocated(neighbours_)) deallocate(neighbours_)
    if (allocated(neighbour_counts_)) deallocate(neighbour_counts_)
  end subroutine freevlist

  !! Calculates the total pair interaction energy of all @p particles. Uses
  !! Verlet neighbour list. 
  !!
  !! @p particles the array of all particles
  !! @p n_particles the number of particles
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
      if(neighbour_counts_(i) == 0) cycle
      do jj = 1, neighbour_counts_(i)
        j = neighbours_(i, jj)
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
  !! @p n_particles the number of particles in @p particles
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
    real(dp) :: distance_moved 
    type(particledat) :: particlevij    
    singleV = 0._dp
    overlap = .false.
    distance_moved = min_distance(simbox, xyz_list_(1:3, i), &
    position(particlei))
    if(distance_moved > 0.5_dp * (r_list_ - cutoff())) then
      call newlist(particles, size(particles), simbox)
    end if
    do j = 1, neighbour_counts_(i)
      particlevij = particles(neighbours_(i, j))
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
  !! @p n_particles the number of particles
  !! @p simbox the simulation cell
  !!
  subroutine newlist(particles, n_particles, simbox)
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    type(poly_box), intent(in) :: simbox
    integer :: i, j
    real(dp) :: r_ij
    do i = 1, n_particles
      neighbour_counts_(i) = 0
      xyz_list_(1, i) = particles(i)%x
      xyz_list_(2, i) = particles(i)%y
      xyz_list_(3, i) = particles(i)%z
    end do 
    do i = 1, n_particles - 1
      do j = i + 1, n_particles
        r_ij = min_distance(simbox, position(particles(i)), &
        position(particles(j)))
        if(r_ij < r_list_ ) then
          neighbour_counts_(i) = neighbour_counts_(i) + 1
          neighbour_counts_(j) = neighbour_counts_(j) + 1 
          neighbours_(i, neighbour_counts_(i)) = j
          neighbours_(j, neighbour_counts_(j)) = i
        end if
      end do
    end do
  end subroutine newlist 

  !! Updates the neighbourlist if needed.
  !!
  !! @p particles the array of particles to be included in the list
  !! @p n_particles the number of particles in @p particles
  !! @p simbox the simulation cell 
  !!
  subroutine updatelist(particles, n_particles, simbox)
    intrinsic max, abs
    type(particledat), dimension(:), intent(in) :: particles     
    integer, intent(in) :: n_particles
    type(poly_box), intent(in) :: simbox
    integer :: i
    do i = 1, n_particles
      if(min_distance(simbox, xyz_list_(1:3, i), position(particles(i))) > &
      0.5_dp * (r_list_ - cutoff())) then
        call newlist(particles, n_particles, simbox)
        exit
      end if
    end do  
  end subroutine updatelist

end module verlet
