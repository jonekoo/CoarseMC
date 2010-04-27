module cell_energy
use cell
use nrtype
use particle
use class_poly_box
use class_pair_potential
use class_nbrcelliterator
use class_parameter_writer
use class_parameterizer
implicit none
private

public :: new_celllist
public :: pairinteractions
public :: new_nbriterator
public :: scaledposition
public :: nbriterator
public :: advance
public :: isdone
public :: value
public :: cell_energy_init
public :: cell_energy_writeparameters
public :: update

real(dp), save :: this_minlength = 7.5_dp
logical, save :: this_iseven = .false.

interface cell_energy_init
  module procedure initwtparameters, initwtreader 
end interface

interface new_celllist
  module procedure new_lists
end interface

interface update
  module procedure cell_energy_update
end interface

interface pairinteractions
  module procedure cell_energy_pairinteractionss, cell_energy_pairinteractionsit
end interface

!! Iterator to iterate through neighbours of a position (a particle).
!! Use new_nbriterator for initialization.
!!
type nbriterator
  private
  type(list), pointer :: cl => NULL()
  type(iterator) :: cit
  type(nbrcelliterator) :: nit
end type

interface advance
  module procedure nbrit_advance
end interface

interface value
  module procedure nbrit_value
end interface

interface isdone
  module procedure nbrit_isdone
end interface

contains

!! Thoughts about keeping the cell list synchronized
!!
!! When do we update the cell list? Whenever a particle or its 
!! neighbours have moved more than 0.5 * (min(cellsides) - rcutoff).
!! Who has this information? The routine that moves the particles. 
!! So that routine could also be responsible of updating the list. 
!! That routine is the sweep routine or the sweep module routines as 
!! before. Any automatic updating would require some sort of bookkeeping
!! in the neighbourlist module or some sort of implementation of the
!! observer pattern. 
!!
!! However, there are many routines that can change the system 
!! configuration. Basicly every change in the system should result in a 
!! test for updating the routine, but the routines that make changes to 
!! the system don't know about each other. They may only know about the 
!! energy calculation routine. If this routine or object is shared among
!! all the updating routines, then it is enough that after each update 
!! to the system, also the energy calculation object is updated. If the 
!! energy calculation routine is not shared, then they can't be 
!! synchronized this way. An observer-observable seems like the only 
!! option. Anyway, the fact that update routines are aware of the need 
!! for update smells like a bad design.
!!
!! Basicly we could also say that every different kind of move has 
!! specific potential energy change related to it. It would be possible
!! to parameterize all the different moves with different potential 
!! energy or "action functional" calculators. 
!!
!! Another, perhaps viable option would be to access particles all the 
!! time through the cell list. In this way it could be automatically 
!! updated and no observing or explicit update mechanism is needed. 
!! Of course this kind of change in the code would send ripples around
!! the system in the form of replacing array index references with some
!! kind of iterators. Here we may come back to the solution that the 
!! particle system should be asked for its energy. Two kind of particle
!! iterators can be used: One which does not give possibility to modify
!! the particles and another which gives that possibility. Actually it 
!! may be that only the iterator for modifying actually needs to be 
!! written. 
!!
!! For optimization of updates, the updating of cell list does not have to 
!! happen when a particle is moved but it may be delayed until the cell list
!! is needed. This may of course require some extra bookkeeping.

!! Constructs a new cell list of @p particles. Given simulation box 
!! dimensions in @p simbox and the minimum cell side length in @p minlength
!! the new list is returned in the variable cl. Note that cell side lengths 
!! may be different in different directions. 
!! 
!! @p cl is the cell list. 
!! @p simbox is the simulation box. 
!! @p particles are the particles to be distributed to the cells.
!! @p nparticles is the number of particles.
!! @p minlength is the minimum length of a cell side. 
!!
function new_lists(simbox, particles) result(cl)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  type(list) :: cl
  integer :: iparticle
  real(dp), dimension(3, size(particles)) :: rs
  integer :: nx, ny, nz
  nx = max(ncells(getx(simbox), this_minlength), 1)
  ny = max(ncells(gety(simbox), this_minlength), 1)
  nz = max(ncells(getz(simbox), this_minlength), 1)
  if (this_iseven) then
    nx = (nx / 2) * 2
    ny = (ny / 2) * 2
    nz = (nz / 2) * 2
  end if 
  do iparticle = 1, size(particles)
    rs(1:3, iparticle) = scaledposition(simbox, particles(iparticle))
  end do
  cl = new_list(rs, nx, ny, nz)
end function

subroutine cell_energy_update(cl, simbox, particles)
  type(list), intent(inout) :: cl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:) :: particles
  cl = new_celllist(simbox, particles)
end subroutine

subroutine initwtreader(parameterreader)
  type(parameterizer), intent(in) :: parameterreader
  call getparameter(parameterreader, 'cellminlength', this_minlength)
  call getparameter(parameterreader, 'isdivisioneven', this_iseven)
end subroutine

subroutine initwtparameters(minlength, iseven)
  real(dp), intent(in) :: minlength
  logical, intent(in) :: iseven
  this_minlength = minlength
  this_iseven = iseven
end subroutine

subroutine cell_energy_writeparameters(writer)
  type(parameter_writer), intent(in) :: writer
  call writecomment(writer, 'cell_energy parameters')
  call writeparameter(writer, 'cellminlength', this_minlength)
  call writeparameter(writer, 'isdivisioneven', this_iseven) 
end subroutine

subroutine cell_energy_pairinteractionss(cl, simbox, particles, i, vpairtot, overlap)
  type(list), intent(in) :: cl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  real(dp), intent(out) :: vpairtot
  logical, intent(out) :: overlap
  type(nbriterator) :: nbrit
  real(dp) :: energy
  nbrit = new_nbriterator(cl, scaledposition(simbox, particles(i)))
  vpairtot = 0._dp
  overlap = .false.
  do 
    if (isdone(nbrit)) exit
    if (value(nbrit) /= i) then
      call pairv(particles(i), particles(value(nbrit)), simbox, energy, &
      overlap)
      if (overlap) return
      vpairtot = vpairtot + energy
    end if
    call advance(nbrit)
  end do
end subroutine

subroutine cell_energy_pairinteractionsit(cl, simbox, particles, vpairtot, overlap)
  type(list), intent(in) :: cl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: vpairtot
  logical, intent(out) :: overlap
  real(dp) :: energy
  integer :: i
  type(nbriterator) :: nbrit
  vpairtot = 0._dp
  overlap = .false.
  do i = 2, size(particles)
    nbrit = new_nbriterator(cl, scaledposition(simbox, particles(i)))
    do 
      if (isdone(nbrit)) exit
      if (value(nbrit) < i) then
        call pairv(particles(i), particles(value(nbrit)), simbox, energy, &
        overlap)
        if (overlap) return
        vpairtot = vpairtot + energy
      end if
      call advance(nbrit)
    end do
  end do
end subroutine

pure function scaledposition(simbox, particle) result(s)
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: particle
  real(dp), dimension(3) :: s
  real(dp), dimension(3) :: r
  r = position(particle)
  s = (/r(1) / getx(simbox), r(2) / gety(simbox), r(3) / getz(simbox)/)
end function

!! Below this nbriterator routines only

!! post-condition: iterator is either isdone or value is valid
!!
function new_nbriterator(cl, scaledpos) result(nbrit)
  type(list), target, intent(in) :: cl
  real(dp), dimension(3), intent(in) :: scaledpos
  type(nbriterator) :: nbrit
  nbrit%cl => cl
  nbrit%nit = new_nbrcelliterator(nx(cl), ny(cl), nz(cl), &
  ix(cl, scaledpos(1)), iy(cl, scaledpos(2)), iz(cl, scaledpos(3)))
  !! :TODO: Think about using the same convention for cell indices in both 
  !! :TODO: interfaces of class_nbrcelliterator and cell.
  nbrit%cit = new_iterator(nbrit%cl, value(nbrit%nit) + 1)
  if (isdone(nbrit%cit)) then
    call advance(nbrit)
  end if
end function

elemental function nbrit_isdone(it) result(done)
  type(nbriterator), intent(in) :: it
  logical :: done
  done = isdone(it%nit)
end function

elemental function nbrit_value(nbrit) result(current)
  type(nbriterator), intent(in) :: nbrit
  integer :: current
  current = value(nbrit%cit)
end function

subroutine nbrit_advance(nbrit) 
  type(nbriterator), intent(inout) :: nbrit
  !! advance in the current cell
  if (.not. isdone(nbrit)) then
    call advance(nbrit%cit)
    do while(isdone(nbrit%cit)) 
      !! advance to next cell
      call advance(nbrit%nit)
      if (isdone(nbrit%nit)) exit
      nbrit%cit = new_iterator(nbrit%cl, value(nbrit%nit) + 1)
    end do
  end if
end subroutine

end module
