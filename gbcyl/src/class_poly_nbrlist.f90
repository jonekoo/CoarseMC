!! This module is a parent module (parent class) for the submodules that 
!! provide neighbourlists. 
module class_poly_nbrlist
  use verlet
  use class_simplelist
  use particle
  use num_kind
  use class_poly_box
  use class_parameter_writer
  use class_parameterizer
  use all_pairs
  use class_pair_potential
  !$ use omp_lib
  implicit none
  private

  public :: pairinteractions
  public :: pnl_writeparameters
  public :: update
  public :: delete
  public :: pnl_init 
  public :: poly_nbrlist
  public :: create_nbrlist
  public :: assignment(=)

  character(len = 80), save :: nbrlisttype = 'stublist'

  type poly_nbrlist
    !private
    type(verletlist), pointer :: vl => NULL()
    type(simplelist), pointer :: sl => NULL()
  end type

  interface pairinteractions
    module procedure pnl_pairinteractions, pnl_totalpairinteractions
  end interface

  interface delete
    module procedure pnl_delete
  end interface

  interface update
    module procedure pnl_update, pnl_updatei
  end interface

  interface assignment(=)
    module procedure assignvl, assignsl, assignpoly
  end interface

  interface nbrmask
    module procedure pnl_nbrmask
  end interface

contains

subroutine pnl_init(parameterreader)
  type(parameterizer), intent(in) :: parameterreader
  call getparameter(parameterreader, 'nbrlisttype', nbrlisttype)
  if (nbrlisttype == 'verletlist') then
    call verlet_init(parameterreader)
  else if (nbrlisttype == 'simplelist') then
    call simplelist_init(parameterreader)
  end if
  call pp_init(parameterreader)
end subroutine

subroutine pnl_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call writeparameter(writer, 'nbrlisttype', nbrlisttype)
  if (nbrlisttype == 'verletlist') call verlet_writeparameters(writer)
  if (nbrlisttype == 'simplelist') call simplelist_writeparameters(writer)
  call pp_writeparameters(writer)
end subroutine

function create_nbrlist(simbox, particles) result(polylist)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  type(poly_nbrlist) :: polylist
  if (nbrlisttype == 'verletlist') then
    polylist = new_verletlist(simbox, particles)
  else if (nbrlisttype == 'simplelist') then
    polylist = new_simplelist(simbox, particles)
  else 
    write(*, *) 'class_poly_nbrlist: Warning! Neighbourlist not created.'
    write(*, *) 'Calculating all pair interactions.'
  end if
end function

subroutine assignvl(polylist, vl)
  type(poly_nbrlist), intent(inout) :: polylist 
  type(verletlist), intent(in) :: vl
  call delete(polylist)
  allocate(polylist%vl)
  polylist%vl = vl
end subroutine

subroutine assignsl(polylist, sl)
  type(poly_nbrlist), intent(inout) :: polylist
  type(simplelist), intent(in) :: sl
  call delete(polylist)
  allocate(polylist%sl)
  polylist%sl=sl
end subroutine

subroutine assignpoly(polylist, another)
  type(poly_nbrlist), intent(inout) :: polylist
  type(poly_nbrlist), intent(in) :: another
  if (associated(another%vl)) then
    polylist = another%vl
  else if (associated(another%sl)) then
    polylist = another%sl
  end if   
end subroutine

subroutine pnl_totalpairinteractions(polylist, simbox, particles, energy, overlap)
  type(poly_nbrlist), intent(in) :: polylist
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap 
  real(dp) :: singleenergy
  logical, dimension(size(particles)) :: mask
  integer :: i
  logical :: overlap_i
  !if (.false.) then
  if (associated(polylist%sl)) then
    call pnl_total_by_cell(polylist, simbox, particles, energy, overlap)
  !$ else if (.true.) then !! Compiled with OpenMP
  !$ overlap = .false.
  !$OMP PARALLEL DO &
  !$OMP DEFAULT(shared) & 
  !$OMP REDUCTION(+:energy) REDUCTION(.or.:overlap) &
  !$OMP private(mask, singleenergy, overlap_i)  
  !$ do i=1, size(particles)
  !$  call nbrmask(polylist, simbox, particles, i, mask)
  !$  call maskedinteractions(mask(i:), simbox, particles(i:), 1, &
  !$       singleenergy, overlap_i)
  !$  overlap = overlap .or. overlap_i
  !$  energy = energy + singleenergy
  !$ end do
  !$OMP END PARALLEL DO
  else
    do i=1, size(particles)
      call nbrmask(polylist, simbox, particles, i, mask)
      call maskedinteractions(mask(i:), simbox, particles(i:), 1, &
      singleenergy, overlap)
      if (overlap) exit
      energy = energy + singleenergy
    end do
  end if
end subroutine


subroutine pnl_total_by_cell(polylist, simbox, particles, &
energy, overlap)
  type(poly_nbrlist), intent(in) :: polylist
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap 
  integer :: i, j, ix, iy, iz
  logical :: mask(size(particles))
  real(dp) :: energy_j
  !integer :: check(size(particles))
  energy = 0._dp
  overlap = .false.
  !! Loop over cells. This can be thought of as looping through a 2 x 2 x 2 
  !! cube
  !! of cells.
  !$OMP PARALLEL default(shared)
  !$OMP DO reduction(+:energy) reduction(.or.:overlap) private(energy_j, i, j, mask) collapse(3)
  do ix=0, polylist%sl%nx-1
  do iy=0, polylist%sl%ny-1
  do iz=0, polylist%sl%nz-1
    !!$ write(*, *) 'thread:', omp_get_thread_num(), 'calculating energy for particles in cell', ix, iy, iz
    call simplelist_nbrmask(polylist%sl, simbox, ix, iy, iz, mask)
    do i=1, polylist%sl%counts(ix, iy, iz)
      j = polylist%sl%indices(i, ix, iy, iz)
      !! There should be a mask by cell -function in simplelist, which could 
      !! be efficiently used in this loop.
      call maskedinteractions(mask(j:), simbox, particles(j:), 1, energy_j, &
      overlap)
      energy = energy + energy_j
      if (overlap) exit
    end do
  end do
  end do
  end do
  !$OMP END DO 
  !$OMP END PARALLEL
end subroutine


pure subroutine maskedinteractions(mask, simbox, particles, i, energy, overlap)
  logical, dimension(:), intent(in) :: mask
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  integer :: j
  real(dp) :: pairenergy
  logical :: overlap_ij
  real(dp) :: rij(3)
  overlap = .false.
  energy = 0._dp
  do j=1, size(particles)
    if(mask(j)) then
      rij = minimage(simbox, position(particles(j)) - position(particles(i)))
      call pairv(particles(i), particles(j), rij, simbox, pairenergy, overlap_ij)
      overlap = overlap .or. overlap_ij
      energy = energy + pairenergy
    end if
  end do
end subroutine

!! Calculate interactions between particles and particle_i
!! 
!! @p particlei must not belong to particles!
subroutine interactions(simbox, particles, particle_i, energy, overlap)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  type(particledat), intent(in) :: particle_i
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  integer :: j
  real(dp) :: pairenergy
  logical :: overlap_ij
  real(dp) :: rij(3)
  overlap = .false.
  energy = 0._dp
  !$OMP PARALLEL
  !$OMP DO PRIVATE(overlap_ij, pairenergy) REDUCTION(+ : energy) REDUCTION(.or. : overlap)
  do j=1, size(particles)
    rij = minimage(simbox, position(particles(j)) - position(particle_i))
    call pairv(particle_i, particles(j), rij, simbox, pairenergy, overlap_ij)
    overlap = overlap .or. overlap_ij
    energy = energy + pairenergy
  end do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine




pure subroutine pnl_pairinteractions(polylist, simbox, particles, i, energy, overlap)
  type(poly_nbrlist), intent(in) :: polylist
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap
  logical :: mask(size(particles))
  call nbrmask(polylist, simbox, particles, i, mask) 
  call maskedinteractions(mask, simbox, particles, i, energy, overlap) 
end subroutine

pure subroutine pnl_nbrmask(polylist, simbox, particles, i, mask) 
  type(poly_nbrlist), intent(in) :: polylist
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  logical, dimension(:), intent(out) :: mask
  if(associated(polylist%vl)) then
    mask = verlet_nbrmask(polylist%vl, size(particles), i)
  else if(associated(polylist%sl)) then
    call simplelist_nbrmask(polylist%sl, simbox, particles, i, mask)
  else 
    mask(:) = .true.
    mask(i) = .false.
  end if
end subroutine

subroutine pnl_update(polylist, simbox, particles)
  type(poly_nbrlist), intent(inout) :: polylist
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  if (associated(polylist%vl)) then
    call update(polylist%vl, simbox, particles)
  else if (associated(polylist%sl)) then
    call update(polylist%sl, simbox, particles)
  end if
end subroutine

subroutine pnl_updatei(polylist, simbox, particles, i)
  type(poly_nbrlist), intent(inout) :: polylist
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  if (associated(polylist%vl)) then
    call update(polylist%vl, simbox, particles, i)
  else if (associated(polylist%sl)) then
    call update(polylist%sl, simbox, particles, i)
  end if
end subroutine


subroutine pnl_delete(polylist)
  type(poly_nbrlist), intent(inout) :: polylist
  if (associated(polylist%vl)) then
    call delete(polylist%vl)
  else if (associated(polylist%sl)) then
    call delete(polylist%sl)
  end if
  nullify(polylist%vl)
  nullify(polylist%sl)
end subroutine

end module
