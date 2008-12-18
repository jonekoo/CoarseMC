!Verlet'n listan käsittelyyn tarvittavia muuttujia ja aliohjelmia. 
module verlet
  use nrtype, only: dp
  use particle, only: particledat, rij, differences, pairV
  implicit none 



  public :: initvlist
  public :: getvlist
  public :: updatelist
  public :: freevlist
  public :: save_state
  public :: load_state
  public :: totpairV



  private

  real(dp), save :: r_list_
  real(dp), save :: r_cutoff_
  integer, save :: n_particles_
  integer, save :: n_neighbours_
  real(dp), dimension(:, :), allocatable, target, save :: xyz_list_
  integer, dimension(:), allocatable, target, save :: neighbour_counts_
  integer, dimension(:, :), allocatable, target, save ::  neighbours_
  integer, save :: n_neighbours_max_
  namelist /verlet_nml/ n_particles_, n_neighbours_, r_list_, r_cutoff_, &
    n_neighbours_max_



  contains

  !Alustaa naapurilistan
  subroutine initvlist(particles, n_particles)
    implicit none
    intrinsic min
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    integer :: astat
    r_list_ = 6.8
    r_cutoff_ = 5.5
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
      write (*, *) 'initvlist: Virhe varattaessa muistia.'
      write (*, *) 'Ohjelman suoritus keskeytyy.'
      stop;
    end if
    call newlist(particles, n_particles)
  end subroutine initvlist



  subroutine freevlist()
    implicit none  
    if (allocated(xyz_list_)) deallocate(xyz_list_)
    if (allocated(neighbours_)) deallocate(neighbours_)
    if (allocated(neighbour_counts_)) deallocate(neighbour_counts_)
  end subroutine freevlist



  subroutine save_state(write_unit)
    implicit none
    integer, intent(in) :: write_unit
    write(write_unit, NML = verlet_nml)
    write(write_unit, *) xyz_list_(1:3, 1:n_particles_)
    write(write_unit, *) neighbour_counts_(1:n_particles_)
    write(write_unit, *) neighbours_(1:n_particles_, 1:n_neighbours_)
  end subroutine save_state



  !! :TODO: probably needs checking if the tables have already been
  !! allocated and to which size they are allocated.
  subroutine load_state(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    read(read_unit, NML = verlet_nml)
    if(.not. allocated(xyz_list_)) allocate(xyz_list_(3, n_particles_))
    read(read_unit, *) xyz_list_(1:3, 1:n_particles_)
    if(.not. allocated(neighbour_counts_)) then 
      allocate(neighbour_counts_(n_particles_))
    end if
    read(read_unit, *) neighbour_counts_(1:n_particles_)
    if(.not. allocated(neighbours_)) then
      allocate(neighbours_(n_particles_, n_neighbours_))
    end if
    read(read_unit, *) neighbours_(1:n_particles_, 1:n_neighbours_)
  end subroutine load_state



  !Palauttaa parivuorovaikutuksen kokonaisenergian
  subroutine totpairV(particles, n_particles, Vtot, ovrlp)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    integer :: i,j,jj
    real(dp), intent(out) :: Vtot 
    real(dp) :: pairE
    logical, intent(out) :: ovrlp
    Vtot = 0.0
    do i=1, n_particles
      if(neighbour_counts_(i) == 0) cycle;
      do jj = 1, neighbour_counts_(i)
        j = neighbours_(i, jj)
        if(j <= i) cycle;
        if (rij(particles(i), particles(j)) < r_cutoff_) then 
          call pairV(particles(i), particles(j), pairE, ovrlp)
          if (ovrlp) return;
          Vtot = Vtot+pairE      
        end if 
      end do
    end do
  end subroutine totpairV



  !! Asettaa osoittimen neighbours_p naapurilistaan, sekä osoittimen 
  !! neighbour_counts_p naapurien lukumäärän sisältävään listaan.
  subroutine getvlist(neighbours_p, neighbour_counts_p)
    implicit none
    integer, dimension(:,:), pointer :: neighbours_p
    integer, dimension(:), pointer :: neighbour_counts_p
    neighbours_p=>neighbours_
    neighbour_counts_p=>neighbour_counts_ 
  end subroutine getvlist



  !Muodostaa uuden naapurilistan hiukkastaulukosta particles
  subroutine newlist(particles, n_particles)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp),dimension(:,:), pointer :: xyz
    integer :: i, j
    xyz=>xyz_list_
    do i = 1, n_particles
      neighbour_counts_(i) = 0
      xyz(1, i) = particles(i)%x
      xyz(2, i) = particles(i)%y
      xyz(3, i) = particles(i)%z
    end do 
    do i = 1, (n_particles-1)
      do j = (i+1), n_particles
         if(rij(particles(i), particles(j)) < r_list_ ) then
           neighbour_counts_(i) = neighbour_counts_(i)+1
           neighbour_counts_(j) = neighbour_counts_(j)+1
           neighbours_(i, neighbour_counts_(i)) = j
           neighbours_(j, neighbour_counts_(j)) = i
         end if
      end do
    end do
  end subroutine newlist 



  !Päivittää naapurilistan tarvittaessa
  subroutine updatelist(particles, n_particles)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles     
    integer, intent(in) :: n_particles
    integer :: i
    real(dp) :: dmax1 = 0.0, dx, dy, dz
    type(particledat) :: oldxyz

    dmax1=0.0
    do i=1, n_particles
      oldxyz%x = xyz_list_(1,i)
      oldxyz%y = xyz_list_(2,i)
      oldxyz%z = xyz_list_(3,i)
      oldxyz%ux = 0.0
      oldxyz%uy = 0.0
      oldxyz%uz = 0.0
      oldxyz%rod = .true.
      call differences(oldxyz, particles(i), dx, dy, dz)
      dx = abs(dx)
      dy = abs(dy)
      dz = abs(dz)
      dmax1 = max(dmax1, dx, dy, dz)
      if(dmax1 > (0.29*(r_list_-r_cutoff_))) then
        call newlist(particles, n_particles)
        exit;
      end if
    end do  
  end subroutine updatelist



end module verlet
