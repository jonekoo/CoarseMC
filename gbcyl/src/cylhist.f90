!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Divides the system in to cylindrical bins and calculates the smctau1, 
!! psi6, smcpsi6 parameters and density for each bin.
!!
!! Usage: The output files for different quantities are given in the file 
!! cylhist.in in the namelist format, e.g.:   
!! &inputnml [binwidth=2.25,] [smctau1file='smctau1.hist',]\
!! [smcpsi6file='smcpsi6.hist',] [densityfile='density.hist',]" \
!! [psi6file='psi6.hist',]/
!!
!! The input configuration(s) are given in the standard input to the program.
!!
!! :TODO: Get rid of magic numbers.
!!
program cylhist
  use state_reader
  use nrtype
  use particle
  use orientational_ordering
  use tau1_module
  use xfunc_module, only: xfunc => rho
  use histogram
  use utils
  use class_factory
  use class_poly_box
  use psi6_module
  use layernormal
  use m_fileunit
  implicit none

  type(particledat), dimension(:), pointer :: particles
  integer :: io_status
  integer :: n_particles
  
  real(dp) :: binwidth = 0.5_dp
  integer :: binningdirection = 1

  !! Cutoff distance for the particles j that are included in the local 
  !! normal of layer around particle i.
  real(dp) :: cutoff = 2._dp 

  !! Only observables which are given a file to write to are 
  !! calculated/written. By default the filenames should be empty
  !! strings.
  character(len=200) :: configurationfile = ''
  character(len=200) :: smcpsi6file = ''
  character(len=200) :: smctau1file = ''
  character(len=200) :: psi6file = ''
  character(len=200) :: layerp2file = ''
  character(len=200) :: ofile = ''
  character(len=200) :: tfile = ''
  character(len=200) :: dfile = ''
  character(len=200) :: nfile = ''

  integer :: configurationunit
  integer :: psi6unit 
  integer :: layerp2unit
  integer :: smctau1unit
  integer :: smcpsi6unit
  integer :: ounit
  integer :: tunit
  integer :: dunit
  integer :: nunit

  namelist /inputnml/ binwidth, cutoff, psi6file, smcpsi6file, smctau1file, &
  & layerp2file, binningdirection, ofile, tfile, dfile, nfile, configurationfile

  integer, dimension(:), pointer :: indices 
  integer :: max_index
  integer :: i_bin
  integer :: allocstat
  integer :: n_bin_particles
  real(dp) :: p2
  real(dp), dimension(3) :: director
  real(sp) :: smctau1
  real(sp) :: tau1
  real(sp) :: layerdistance

  !! Parameters to control which observables are calculated and written.
  logical :: writepsi6 = .false.
  logical :: writesmcpsi6 = .false.
  logical :: writesmctau1 = .false.
  logical :: writelayerp2 = .false.
  logical :: writeorientation = .false.
  logical :: writedensity = .false.
  logical :: writetau1 = .false.
  logical :: write_particle_count = .false.

  type(poly_box) :: simbox
  type(factory) :: afactory
  integer, parameter :: stdin = 5
 
  !! Read input parameters
  read(stdin, NML = inputnml) 

  !! Open configurationfile
  if (trim(configurationfile) /= '') then
    configurationunit=fileunit_getfreeunit()
    open(unit=configurationunit, file=configurationfile, action='READ', status='OLD')
  else 
    write(*, *) 'Set configurationfile!' 
    stop
  end if

  if (trim(psi6file) /= '') then
    call open_write_unit(psi6file, psi6unit)
    writepsi6 = .true.
  end if
  if (trim(layerp2file) /= '') then
    call open_write_unit(layerp2file, layerp2unit)
    writelayerp2 = .true.
  end if
  if (trim(smctau1file) /= '') then
    call open_write_unit(smctau1file, smctau1unit)
    writesmctau1 = .true.
  end if
  if (trim(smcpsi6file) /= '') then
    call open_write_unit(smcpsi6file, smcpsi6unit)
    writesmcpsi6 = .true.
  end if
  if (trim(ofile) /= '') then
    call open_write_unit(ofile, ounit)
    writeorientation = .true.
  end if
  if (trim(tfile) /= '') then
    call open_write_unit(tfile, tunit)
    writetau1 = .true.
  end if
  if (trim(dfile) /= '') then
    call open_write_unit(dfile, dunit)
    writedensity = .true.
  end if
  if (trim(nfile) /= '') then
    call open_write_unit(nfile, nunit)
    write_particle_count = .true.
  else 
    stop 'file for particle number not given'
  end if

  do  
    call readstate(afactory, configurationunit, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    end if
    n_particles = size(particles)
    if (associated(indices)) then
      if(size(indices) /= n_particles) then
        if(associated(indices)) deallocate(indices)
        allocate(indices(n_particles), stat=allocstat)
        if(allocstat /= 0) stop 'Allocation of indices failed.'
      end if
    else 
      allocate(indices(n_particles), stat=allocstat)
      if(allocstat /= 0) stop 'Allocation of indices failed.'
    end if

    !! Solve bin indices for all the particles
    if (binningdirection == -1) then
      call bin_indices(particles, n_particles, xfunc, binwidth, indices,&
      getx(simbox)/2._dp, binningdirection)
    else 
      call bin_indices(particles, n_particles, xfunc, binwidth, indices)
    end if
    max_index = maxval(indices)
    do i_bin = 1, max_index
      n_bin_particles = count(indices == i_bin)

      write(nunit, '('//fmt_char_int()//',1X)', ADVANCE='NO') n_bin_particles

      if (writepsi6) then
        if(n_bin_particles == 0) then         
          write(psi6unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write psi6 profile
          call orientation_parameter(pack(particles, (indices == i_bin)), &
          n_bin_particles, p2, director) 
          write(psi6unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          psi6_bulk_masked(simbox, particles, (indices == i_bin), director) 
        end if
      end if

      if (writesmcpsi6 .or. writesmctau1 .or. writelayerp2) then
        if(n_bin_particles == 0) then        
          write(smcpsi6unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate layer normal for i_bin
          call globalnormal(simbox, pack(particles(1:n_particles), &
          indices == i_bin), cutoff, p2, director)
        end if
      end if

      if (writelayerp2) then
        if(n_bin_particles == 0) then    
          write(layerp2unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Write the orientation parameter for layer normals.
          write(layerp2unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') p2
        end if
      end if

      if (writesmcpsi6) then
        if(n_bin_particles == 0) then    
          write(smcpsi6unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write psi6 profile with respect to averaged layer 
          !! normal for i_bin.
          write(smcpsi6unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          psi6_bulk_masked(simbox, particles, (indices == i_bin), director) 
        end if
      end if

      if (writesmctau1) then
        if(n_bin_particles == 0) then    
          write(smctau1unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write smctau1 profile        
          call tau1_routine(pack(particles(1:n_particles), indices == i_bin), &
          real(director, sp), smctau1, layerdistance)
          write(smctau1unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') smctau1
        end if
      end if
    
      if (writedensity) then
        !! Calculate and write density profile
        write(dunit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
        real(n_bin_particles, dp)/binvolume(i_bin, getx(simbox)/2._dp, &
        getz(simbox), binningdirection)
      end if

      if (writeorientation) then
        if(n_bin_particles == 0) then    
          write(ounit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write orientation parameter profile
          call orientation_parameter(pack(particles(1:n_particles), &
          indices == i_bin), n_bin_particles, p2, director)
          write(ounit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') p2 
        end if
      end if

      if (writetau1) then
        if(n_bin_particles == 0) then    
          write(tunit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write tau1 profile        
          call tau1_routine(pack(particles(1:n_particles), indices == i_bin), &
          real(director, sp), tau1, layerdistance)
          write(tunit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') tau1
        end if
      end if

    end do

    !! Append a newline after each row of bins
    if (writepsi6) then
      write(psi6unit, *) ''
    end if
    if (writesmctau1) then
      write(smctau1unit, *) ''
    end if
    if (writesmcpsi6) then
      write(smcpsi6unit, *) ''
    end if
    if (writelayerp2) then
      write(layerp2unit, *) ''
    end if
    !! Append a newline after each row of bins
    if (writeorientation) then
      write(ounit, *) ''
    end if
    if (writetau1) then
      write(tunit, *) ''
    end if
    if (writedensity) then
      write(dunit, *) ''
    end if
    write(nunit, *) ''
  end do

  call finalize

  contains

  function binvolume(i_bin, maxr, height, binningdirection)
    integer, intent(in) :: i_bin
    real(dp), intent(in) :: maxr
    real(dp), intent(in) :: height
    integer, intent(in) :: binningdirection
    real(dp) :: binvolume
    if (binningdirection == -1) then
      binvolume = 4._dp * atan(1._dp) * height * ((maxr - real(i_bin-1, dp) * &
      binwidth)**2 - max(maxr - real(i_bin, dp) * binwidth, 0._dp)**2)    
    else
      binvolume = 4._dp * atan(1._dp) * height * ((max(real(i_bin, dp) * &
      binwidth, maxr))**2 - (real(i_bin-1, dp) * binwidth)**2)    
    end if    
  end function 

  subroutine finalize
    if(associated(particles)) deallocate(particles)
    if(associated(indices)) deallocate(indices)
    close(configurationunit)
    if (writepsi6) then
      close(psi6unit)
    end if
    if (writesmctau1) then
      close(smctau1unit)
    end if
    if (writesmcpsi6) then
      close(smcpsi6unit)
    end if
    if (writelayerp2) then
      close(layerp2unit)
    end if
  end subroutine

  subroutine open_write_unit(file_name, write_unit)
    implicit none
    character(len = *), intent(in) :: file_name
    integer, intent(out) :: write_unit
    integer :: opened
    write_unit=fileunit_getfreeunit()
    open(write_unit, FILE = file_name, status = 'replace', &
      position = 'append', iostat = opened)
    if(opened/=0) then
      write(*,*) 'Opening of file ', file_name, ' failed.'
      stop 
    end if
  end subroutine open_write_unit
end program

