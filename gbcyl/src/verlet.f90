module verlet
  !! Module for calculating pair interactions of particles using a Verlet 
  !! neighbour list. 
  use nrtype, only: dp
  use particle, only: particledat, position, getmaxmoves
  use class_poly_box
  use class_parameterizer
  use class_parameter_writer
  implicit none 
  private  

  public :: verlet_init
  public :: verlet_writeparameters
  public :: new_verletlist
  public :: verletlist
  public :: update
  public :: delete
  public :: nbrmask
 
  real(dp), save :: rlist = 6.8_dp

  !! When any cartesian coordinate of any particle has changed more than
  real(dp), save :: updatethreshold = 0.377_dp
  !! the neighbourlist will be updated by verlet_updatei. 

  integer, save :: nneighboursmax = 500

  type verletlist
    private
    integer :: nneighbours = 500
    real(dp), dimension(:, :), pointer :: xyzlist => NULL()
    integer, dimension(:), pointer :: neighbourcounts => NULL()
    integer, dimension(:, :), pointer ::  neighbours => NULL()
  end type 

  interface verlet_init
    module procedure verlet_initwtparameterizer, verlet_initwtvalues
  end interface

  interface delete
    module procedure verlet_delete
  end interface

  interface update
    module procedure verlet_update, verlet_updatei
  end interface

  interface nbrmask
    module procedure verlet_nbrmask
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
    call getparameter(reader, 'n_neighbours_max', nneighboursmax)
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


  subroutine verlet_updatei(vl, simbox, particles, i)
    type(verletlist), intent(inout) :: vl
    type(poly_box), intent(in) :: simbox
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: i 
    if (any(minimage(simbox, position(particles(i))-vl%xyzlist(:,i)) > updatethreshold)) then
      call update(vl, simbox, particles)
    end if
  end subroutine


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
        rijvec = minimage(simbox, position(particles(j)) - &
        position(particles(i)))
        if(dot_product(rijvec, rijvec) < rlist*rlist ) then
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

  pure function verlet_nbrmask(vl, nparticles, i)
    type(verletlist), intent(in) :: vl
    integer, intent(in) :: nparticles
    integer, intent(in) :: i
    logical, dimension(nparticles) :: verlet_nbrmask
    integer :: j
    verlet_nbrmask(:) = .false.
    do j = 1, vl%neighbourcounts(i)
      verlet_nbrmask(vl%neighbours(i, j)) = .true.
    end do
  end function

end module verlet
