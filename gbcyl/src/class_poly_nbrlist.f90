!! This module is a parent module (parent class) for the submodules that 
!! provide neighbourlists. 
module class_poly_nbrlist
  use verlet
  !use cell_energy
  !use cell
  use class_simplelist
  use particle
  use nrtype
  use class_poly_box
  use class_parameter_writer
  use class_parameterizer
  use all_pairs
  use class_pair_potential
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
    private
    type(verletlist), pointer :: vl => NULL()
    !type(list), pointer :: cl => NULL()
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
!  if (nbrlisttype == 'celllist') then
!    call cell_energy_init(parameterreader)
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
!  if (nbrlisttype == 'celllist') call cell_energy_writeparameters(writer)
  if (nbrlisttype == 'simplelist') call simplelist_writeparameters(writer)
  call pp_writeparameters(writer)
end subroutine

function create_nbrlist(simbox, particles) result(polylist)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  type(poly_nbrlist) :: polylist
  if (nbrlisttype == 'verletlist') then
    polylist = new_verletlist(simbox, particles)
!  else if (nbrlisttype == 'celllist') then
!    polylist = new_celllist(simbox, particles)
  else if (nbrlisttype == 'simplelist') then
    polylist = new_simplelist(simbox, particles)
  else 
    write(*, *) 'class_poly_nbrlist: Warning! Neighbourlist not created.'
    write(*, *) 'Calculating all pair interactions.'
  end if
end function

!subroutine assigncl(polylist, cl)
!  type(poly_nbrlist), intent(inout) :: polylist 
!  type(list), intent(in) :: cl
!  call delete(polylist)
!  allocate(polylist%cl)
!  polylist%cl = cl
!end subroutine

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
  !else if (associated(another%cl)) then
  !  polylist = another%cl
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
  do i=1, size(particles)
    mask = nbrmask(polylist, simbox, particles, i)
    call maskedinteractions(mask(i:), simbox, particles(i:), 1, &
      singleenergy, overlap)
    if (overlap) exit
    energy = energy + singleenergy
  end do
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
  energy = 0._dp
  !$OMP DO DEFAULT(shared) PRIVATE(j, pairenergy), REDUCTION(+:energy, .OR. overlap)
  do j=1, size(particles)
    if(mask(j)) then
      call pairv(particles(i), particles(j), simbox, pairenergy, overlap)
      if (overlap) exit
      energy = energy + pairenergy
    end if
  end do
  !$OMP END DO 
end subroutine

subroutine pnl_pairinteractions(polylist, simbox, particles, i, energy, overlap)
  type(poly_nbrlist), intent(in) :: polylist
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  real(dp), intent(out) :: energy
  logical, intent(out) :: overlap 
  call maskedinteractions(nbrmask(polylist, simbox, particles, i), simbox, particles, i, energy, overlap) 
end subroutine

function pnl_nbrmask(polylist, simbox, particles, i) result(mask)
  type(poly_nbrlist), intent(in) :: polylist
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  logical, dimension(size(particles)) :: mask
!  if (associated(polylist%cl)) then
!    mask = nbrmask(polylist%cl, simbox, particles, i)
  if(associated(polylist%vl)) then
    mask = nbrmask(polylist%vl, size(particles), i)
  else if(associated(polylist%sl)) then
    mask = nbrmask(polylist%sl, simbox, particles, i)
  else 
    mask(:) = .true.
    mask(i) = .false.
  end if
end function

subroutine pnl_update(polylist, simbox, particles)
  type(poly_nbrlist), intent(inout) :: polylist
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
!  if (associated(polylist%cl)) then
!    call update(polylist%cl, simbox, particles)
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
  !if (associated(polylist%cl)) then
  !  call update(polylist%cl, simbox, particles, i)
  if (associated(polylist%vl)) then
    call update(polylist%vl, simbox, particles, i)
  else if (associated(polylist%sl)) then
    call update(polylist%sl, simbox, particles, i)
  end if
end subroutine


subroutine pnl_delete(polylist)
  type(poly_nbrlist), intent(inout) :: polylist
  !if (associated(polylist%cl)) then
  !  call delete(polylist%cl)
  if (associated(polylist%vl)) then
    call delete(polylist%vl)
  else if (associated(polylist%sl)) then
    call delete(polylist%sl)
  end if
  !nullify(polylist%cl)
  nullify(polylist%vl)
  nullify(polylist%sl)
end subroutine

end module
