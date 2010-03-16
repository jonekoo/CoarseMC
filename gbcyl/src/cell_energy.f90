module cell_energy
use cell
use nrtype
use particle
use class_poly_box
use class_pair_potential
implicit none
private

public :: new_list
public :: pair_interactions
public :: new_nbr_iterator
public :: scaled_position
public :: nbr_iterator

interface new_list
  module procedure new_list_s
end interface

interface pair_interactions
  module procedure pair_interactions_s, pair_interactions_it
end interface

!! Iterator to iterate through neighbours of a position (a particle).
!! Use new_nbr_iterator for initialization.
!!
type nbr_iterator
  private
  type(list), pointer :: cl => NULL()
  type(iterator) :: cit
  integer :: ix, iy, iz
  integer :: dix, diy, diz
  logical :: is_done = .true. 
end type

interface advance
  module procedure nbrit_advance
end interface

interface value
  module procedure nbrit_value
end interface

interface is_done
  module procedure nbrit_is_done
end interface

contains

!! Thoughts about keeping the cell list synchronized
!!
!! When do we update the cell list? Whenever a particle or its 
!! neighbours have moved more than 0.5 * (min(cell_sides) - r_cutoff).
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
!! dimensions in @p simbox and the minimum cell side length in @p min_length
!! the new list is returned in the variable cl. Note that cell side lengths 
!! may be different in different directions. 
!! 
!! @p cl is the cell list. 
!! @p simbox is the simulation box. 
!! @p particles are the particles to be distributed to the cells.
!! @p n_particles is the number of particles.
!! @p min_length is the minimum length of a cell side. 
!!
pure function new_list_s(simbox, particles, min_length) result(cl)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(in) :: min_length
  type(list) :: cl
  integer :: i_particle
  real(dp), dimension(3, size(particles)) :: rs
  integer :: nx, ny, nz
  nx = n_cells(get_x(simbox), min_length)
  ny = n_cells(get_y(simbox), min_length)
  nz = n_cells(get_z(simbox), min_length)
  do i_particle = 1, size(particles)
    rs(1:3, i_particle) = scaled_position(simbox, particles(i_particle))
  end do
  cl = new_list_f(rs, nx, ny, nz)
end function

subroutine pair_interactions_s(cl, simbox, particles, i, vpairtot, overlap)
  type(list), intent(in) :: cl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  integer, intent(in) :: i
  real(dp), intent(out) :: vpairtot
  logical, intent(out) :: overlap
  type(nbr_iterator) :: nbrit
  real(dp) :: energy
  nbrit = new_nbr_iterator(cl, scaled_position(simbox, particles(i)))
  vpairtot = 0._dp
  overlap = .false.
  do 
    if (is_done(nbrit)) exit
    if (value(nbrit) /= i) then
      call pairv(particles(i), particles(value(nbrit)), simbox, energy, &
      overlap)
      if (overlap) return
      vpairtot = vpairtot + energy
    end if
  end do
end subroutine

subroutine pair_interactions_it(cl, simbox, particles, vpairtot, overlap)
  type(list), intent(in) :: cl
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: vpairtot
  logical, intent(out) :: overlap
  real(dp) :: energy
  integer :: i
  type(nbr_iterator) :: nbrit
  vpairtot = 0._dp
  overlap = .false.
  do i = 2, size(particles)
    nbrit = new_nbr_iterator(cl, scaled_position(simbox, particles(i)))
    do 
      if (is_done(nbrit)) exit
      if (value(nbrit) < i) then
        call pairv(particles(i), particles(value(nbrit)), simbox, energy, &
        overlap)
        if (overlap) return
        vpairtot = vpairtot + energy
      end if
    end do
  end do
end subroutine


!! Don't use this!
!!
subroutine pair_interactions_t(simbox, particles, vpairtot, overlap)
  type(poly_box), intent(in) :: simbox
  type(particledat), dimension(:), intent(in) :: particles
  real(dp), intent(out) :: vpairtot
  logical, intent(out) :: overlap
  type(list) :: clist
  integer :: i_cell
  type(iterator) :: it1
  type(iterator) :: it2
  real(dp) :: pot_energy
  integer :: ix1
  integer :: iy1
  integer :: iz1
  integer :: ix2
  integer :: iy2
  integer :: iz2
  integer :: dix2
  integer :: diy2
  integer :: diz2
  vpairtot = 0._dp
  overlap = .false.
  if (size(particles) == 0) return
  clist = new_list(simbox, particles, 2._dp * max_trans(particles(1)) + cutoff())
  do ix1 = 0, nx_(clist) - 1
  do iy1 = 0, ny_(clist) - 1
  do iz1 = 0, nz_(clist) - 1
    it1 = new_iterator(clist, i_cell)
    do 
      if (is_done(it1)) exit
      do dix2 = -1, 1
      do diy2 = -1, 1
      do diz2 = -1, 1
        ix2 = ix1 + dix2 
        iy2 = iy1 + diy2
        iz2 = iz2 + diz2
        if (is_xperiodic(simbox)) ix2 = mod(ix2 + nx_(clist), nx_(clist))
        if (is_yperiodic(simbox)) iy2 = mod(iy2 + ny_(clist), ny_(clist))
        if (is_zperiodic(simbox)) iz2 = mod(iz2 + nz_(clist), nz_(clist))
        it2 = new_iterator(clist, cell_index(clist, ix2, iy2, iz2))
        do 
          if (is_done(it2)) exit
          if (value(it1) < value(it2)) then
            call pairv(particles(value(it1)), particles(value(it2)), simbox, pot_energy, overlap)
            if (overlap) return
            vpairtot = vpairtot + pot_energy
          end if
          call advance(it2)
        end do
      end do
      end do
      end do
      call advance(it1)
    end do
  end do
  end do
  end do 
end subroutine

function new_nbr_iterator(cl, position) result(nit)
  type(list), target, intent(in) :: cl
  real(dp), dimension(3), intent(in) :: position
  type(nbr_iterator) :: nit
  nit%cl => cl
  nit%ix = ix(cl, position(1))
  nit%iy = iy(cl, position(2))
  nit%iz = iz(cl, position(3))
  nit%dix = -1
  nit%diy = -1
  nit%diz = -1
  nit%is_done = .false.
end function

pure subroutine nbrit_advance(nit)
  type(nbr_iterator), intent(inout) :: nit
  if (is_done(nit)) then 
    return
  else if (is_done(nit%cit)) then
    if (nit%dix < 1) then 
      nit%dix = nit%dix + 1
    else if (nit%diy < 1) then 
      nit%diy = nit%diy + 1
      nit%dix = -1
    else if (nit%diz < 1) then
      nit%diz = nit%diz + 1
      nit%dix = -1
      nit%diy = -1
    else
      nit%is_done = .true.
      return
    end if
    nit%cit = new_iterator(nit%cl, cell_index(nit%cl, nit%ix + nit%dix, nit%iy + nit%diy, nit%iz + nit%diz)) 
  else
    call advance(nit%cit)
  end if 
end subroutine

elemental function nbrit_is_done(nit) result(done)
  type(nbr_iterator), intent(in) :: nit
  logical :: done
  done = nit%is_done
end function

elemental function nbrit_value(nit) result(current)
  type(nbr_iterator), intent(in) :: nit
  integer :: current
  current = value(nit%cit)
end function

pure function scaled_position(simbox, particle) result(s)
  type(poly_box), intent(in) :: simbox
  type(particledat), intent(in) :: particle
  real(dp), dimension(3) :: s
  real(dp), dimension(3) :: r
  r = position(particle)
  s = (/r(1)/get_x(simbox), r(2)/get_y(simbox), r(3)/get_z(simbox)/)
end function

end module
