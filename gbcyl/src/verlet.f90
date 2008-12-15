!Verlet'n listan käsittelyyn tarvittavia muuttujia ja aliohjelmia. 
module verlet
  use nrtype, only: dp
  use particle, only: particledat, rij, differences
  implicit none 

  real(dp), parameter :: rlist = 6.8, rcut = 5.5
  integer, private, save :: n_particles_
  integer, private, save :: n_neighbours_
  real(dp), dimension(:, :), allocatable, target, private, save :: xyzlist
  integer, dimension(:), allocatable, target, private, save :: Nvlist
  integer, dimension(:, :), allocatable, target, private, save ::  vlist
  namelist /verlet_nml/ n_particles_, n_neighbours_
  PRIVATE :: verlet_nml



  contains



  !Alustaa naapurilistan
  subroutine initvlist(particles, n_particles)
    implicit none
    intrinsic min
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    integer :: astat
    integer, parameter :: n_neighbours_max = 500
    n_particles_ = n_particles
    n_neighbours_ = min(n_particles_, n_neighbours_max) 
    if (n_particles_ < n_neighbours_) then
      n_neighbours_ = n_particles_
    end if
    allocate(xyzlist(3, n_particles_), vlist(n_particles_, n_neighbours_), & 
           & Nvlist(n_particles_), stat = astat)
    if (astat /= 0) then
      write (*, *) 'initvlist: Virhe varattaessa muistia.'
      write (*, *) 'Ohjelman suoritus keskeytyy.'
      stop;
    end if
    call newlist(particles, n_particles)
  end subroutine initvlist



  subroutine freevlist()
    implicit none  
    if (allocated(xyzlist)) deallocate(xyzlist)
    if (allocated(vlist)) deallocate(vlist)
    if (allocated(Nvlist)) deallocate(Nvlist)
  end subroutine freevlist



  subroutine save_state(write_unit)
    implicit none
    integer, intent(in) :: write_unit
    write(write_unit, NML = verlet_nml)
    write(write_unit, *) xyzlist(1:3, 1:n_particles_)
    write(write_unit, *) Nvlist(1:n_particles_)
    write(write_unit, *) vlist(1:n_particles_, 1:n_neighbours_)
  end subroutine save_state



  !! :TODO: probably needs checking if the tables have already been
  !! allocated and to which size they are allocated.
  subroutine load_state(read_unit)
    implicit none
    integer, intent(in) :: read_unit
    read(read_unit, NML = verlet_nml)
    if(.not. allocated(xyzlist)) allocate(xyzlist(3, n_particles_))
    read(read_unit, *) xyzlist(1:3, 1:n_particles_)
    if(.not. allocated(Nvlist)) allocate(Nvlist(n_particles_))
    read(read_unit, *) Nvlist(1:n_particles_)
    if(.not. allocated(vlist)) allocate(vlist(n_particles_, n_neighbours_))
    read(read_unit, *) vlist(1:n_particles_, 1:n_neighbours_)
  end subroutine load_state



  !Asettaa osoittimen vlistp naapurilistaan, sekä osoittimen Nvlistp
  !naapurien lukumäärän sisältävään listan
  subroutine getvlist(vlistp, Nvlistp)
    implicit none
    integer, dimension(:,:), pointer :: vlistp
    integer, dimension(:), pointer :: Nvlistp
    vlistp=>vlist
    Nvlistp=>Nvlist 
  end subroutine getvlist



  !Muodostaa uuden naapurilistan hiukkastaulukosta particles
  subroutine newlist(particles, n_particles)
    implicit none
    type(particledat), dimension(:), intent(in) :: particles
    integer, intent(in) :: n_particles
    real(dp),dimension(:,:), pointer :: xyz
    integer :: i, j
    xyz=>xyzlist
    do i = 1, n_particles
      Nvlist(i) = 0
      xyz(1, i) = particles(i)%x
      xyz(2, i) = particles(i)%y
      xyz(3, i) = particles(i)%z
    end do 
    do i = 1, (n_particles-1)
      do j = (i+1), n_particles
         if(rij(particles(i), particles(j)) < rlist ) then
           Nvlist(i) = Nvlist(i)+1
           Nvlist(j) = Nvlist(j)+1
           vlist(i, Nvlist(i)) = j
           vlist(j, Nvlist(j)) = i
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
      oldxyz%x = xyzlist(1,i)
      oldxyz%y = xyzlist(2,i)
      oldxyz%z = xyzlist(3,i)
      oldxyz%ux = 0.0
      oldxyz%uy = 0.0
      oldxyz%uz = 0.0
      oldxyz%rod = .true.
      call differences(oldxyz, particles(i), dx, dy, dz)
      dx = abs(dx)
      dy = abs(dy)
      dz = abs(dz)
      dmax1 = max(dmax1, dx, dy, dz)
      if(dmax1 > (0.29*(rlist-rcut))) then
        call newlist(particles, n_particles)
        exit;
      end if
    end do  
  end subroutine updatelist



end module verlet
