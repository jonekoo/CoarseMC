!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Divides the system in to cylindrical bins and calculates for each bin
!! - Number of particles
!! - Number density of particles
!!and 
!! - Orientation parameter of particles with respect to cylinder axis
!! - Translational order with respect to the layering of the whole system tau1
!! - Short range hexagonal order of the neighbouring particles psi6
!! For all those above the reference direction is the cylinder axis.
!!
!! To recognize smectic C features one can also calculate
!! - The orientation parameter of the smectic layer layerp2
!! - Translational order smctau1 and 
!! - short range hexagonal order smcpsi6. 
!! These latter ones have the average layer normal as their reference direction.
!! 
!! Usage: The output files for different quantities are given in standard input
!! in the F95 namelist format as in the example below:
!!
!! &inputnml configuration_file='configurations.0', n_file, bin_width=0.25, \
!! n_bins=36, [binning_direction=-1,] \
!! [smctau1_file='smctau1.hist',] [smcpsi6_file='smcpsi6.hist',] \
!! [density_file='density.hist',] [psi6_file='psi6.hist',] \
!! [p2_file='p2.hist',] [tau1_file='tau1.hist',]/
!!
!! A parameter is calculated and printed only if its output file is given. 
!! binning_direction=-1 means the division to bins is started from the 
!! cylinder wall.
!! 
program cylhist
  use state_reader
  use nrtype
  use particle
  use orientational_ordering
  use tau1_module
  use tau1_negative
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
  
  real(dp) :: bin_width = 0.5_dp
  integer :: binning_direction = 1

  !! Cutoff distance for the particles j that are included in the local 
  !! normal of layer around particle i.
  real(dp) :: cutoff = 2._dp 

  !! Only observables which are given a file to write to are 
  !! calculated/written. By default the filenames should be empty
  !! strings.
  character(len=200) :: configuration_file = ''
  character(len=200) :: smcpsi6_file = ''
  character(len=200) :: smctau1_file = ''
  character(len=200) :: psi6_file = ''
  character(len=200) :: layerp2_file = ''
  character(len=200) :: p2_file = ''
  character(len=200) :: tau1_file = ''
  character(len=200) :: density_file = ''
  character(len=200) :: n_file = ''

  integer :: configuration_unit
  integer :: psi6_unit 
  integer :: layerp2_unit
  integer :: smctau1_unit
  integer :: smcpsi6_unit
  integer :: p2_unit
  integer :: tau1_unit
  integer :: density_unit
  integer :: n_unit

  namelist /inputnml/ bin_width, cutoff, psi6_file, smcpsi6_file, smctau1_file, &
  & layerp2_file, binning_direction, p2_file, tau1_file, density_file, n_file, configuration_file, n_bins

  integer, dimension(:), pointer :: indices 
  integer :: n_bins
  integer :: i_bin
  integer :: allocstat
  integer :: n_bin_particles
  real(dp) :: tensor(3,3)
  real(dp) :: p2
  real(dp) :: layer_p2
  real(dp) :: layer_normal(3)
  real(sp) :: bin_layer_normal(3) 
  real(sp) :: tau1_value
  real(sp) :: smctau1_value
  real(sp) :: layer_distance
  real(sp) :: smclayer_distance
  real(sp), dimension(3), parameter :: z_axis=(/0._sp, 0._sp, 1._sp/)

  !! Parameters to control which observables are calculated and written.
  logical :: write_psi6 = .false.
  logical :: write_smcpsi6 = .false.
  logical :: write_smctau1 = .false.
  logical :: write_layerp2 = .false.
  logical :: write_p2 = .false.
  logical :: write_density = .false.
  logical :: write_tau1 = .false.
  logical :: write_particle_count = .false.

  type(poly_box) :: simbox
  type(factory) :: afactory
  integer, parameter :: stdin = 5
 
  !! Read input parameters
  read(stdin, NML = inputnml) 

  if (n_bins==0) stop 'Give n_bins!'

  !! Open configuration_file
  if (trim(configuration_file) /= '') then
    configuration_unit=fileunit_getfreeunit()
    open(unit=configuration_unit, file=configuration_file, action='READ', status='OLD')
  else 
    write(*, *) 'Set configuration_file!' 
    stop
  end if

  if (trim(psi6_file) /= '') then
    call open_write_unit(psi6_file, psi6_unit)
    write_psi6 = .true.
  end if
  if (trim(layerp2_file) /= '') then
    call open_write_unit(layerp2_file, layerp2_unit)
    write_layerp2 = .true.
  end if
  if (trim(smctau1_file) /= '') then
    call open_write_unit(smctau1_file, smctau1_unit)
    write_smctau1 = .true.
  end if
  if (trim(smcpsi6_file) /= '') then
    call open_write_unit(smcpsi6_file, smcpsi6_unit)
    write_smcpsi6 = .true.
  end if
  if (trim(p2_file) /= '') then
    call open_write_unit(p2_file, p2_unit)
    write_p2 = .true.
  end if
  if (trim(tau1_file) /= '') then
    call open_write_unit(tau1_file, tau1_unit)
    write_tau1 = .true.
  end if
  if (trim(density_file) /= '') then
    call open_write_unit(density_file, density_unit)
    write_density = .true.
  end if
  if (trim(n_file) /= '') then
    call open_write_unit(n_file, n_unit)
    write_particle_count = .true.
  else 
    stop 'n_file not given'
  end if

  do  
    call readstate(afactory, configuration_unit, simbox, particles, io_status)
    if (io_status < 0) then
      exit
    end if
    n_particles = size(particles)

    if (write_p2) then
      !! Calculate and save orientation parameter for the whole system with 
      !! respect to z-axis.
      call orientation_tensor(particles, n_particles, tensor)
      write(p2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') tensor(3,3) 
    end if

    if (write_tau1) then
      !! Calculate and write global tau1
      call tau1_routine(particles, z_axis, tau1_value, layer_distance)
      write(tau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') tau1_value
    end if

    if (write_psi6) then
      write(psi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
      psi6_bulk(simbox, particles, real(z_axis, dp)) 
    end if

    if (write_smcpsi6 .or. write_smctau1 .or. write_layerp2) then
      !! Calculate average orientation of layer normals
      call globalnormal(simbox, particles, cutoff, layer_p2, layer_normal)
    end if

    if (write_layerp2) then
      write(layerp2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') layer_p2
    end if

    if (write_smctau1) then
      !! Calculate and write global tau1
      call tau1_routine(particles, real(layer_normal, sp), tau1_value, smclayer_distance)
      write(smctau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') smctau1_value
    end if

    if (write_smcpsi6) then
      write(smcpsi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
      psi6_bulk(simbox, particles, real(layer_normal, dp)) 
    end if

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
    if (binning_direction == -1) then
      call bin_indices(particles, n_particles, xfunc, bin_width, indices,&
      getx(simbox)/2._dp, binning_direction)
    else 
      call bin_indices(particles, n_particles, xfunc, bin_width, indices)
    end if
    do i_bin = 1, n_bins
      n_bin_particles = count(indices == i_bin)

      write(n_unit, '('//fmt_char_int()//',1X)', ADVANCE='NO') n_bin_particles

      if (write_density) then
        !! Calculate and write density profile
        write(density_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
        real(n_bin_particles, dp)/binvolume(i_bin, getx(simbox)/2._dp, &
        getz(simbox), binning_direction)
      end if

      if (write_p2) then
        if(n_bin_particles == 0) then    
          write(p2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write orientation parameter profile
          call orientation_tensor(pack(particles(1:n_particles), &
          & indices==i_bin), n_bin_particles, tensor)
          write(p2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') tensor(3,3) 
        end if
      end if

      if (write_tau1) then
        if(n_bin_particles == 0) then    
          write(tau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write tau1 profile        
          write(tau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          tau1(pack(particles(1:n_particles), indices==i_bin), z_axis, layer_distance)
        end if
      end if

      if (write_psi6) then
        if(n_bin_particles == 0) then         
          write(psi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write psi6 profile
          write(psi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          psi6_bulk_masked(simbox, particles, (indices == i_bin), real(z_axis, dp)) 
        end if
      end if

      if (write_layerp2) then
        if(n_bin_particles == 0) then    
          write(layerp2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Write the orientation parameter for layer normals.
          write(layerp2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          dot_product(bin_layer_normal, layer_normal)
        end if
      end if

      if (write_smctau1) then
        if(n_bin_particles == 0) then    
          write(smctau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write smctau1 profile        
          write(smctau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          tau1(pack(particles(1:n_particles), indices==i_bin), real(layer_normal, sp), smclayer_distance)
        end if
      end if
    
      if (write_smcpsi6) then
        if(n_bin_particles == 0) then    
          write(smcpsi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 0._dp
        else
          !! Calculate and write psi6 profile with respect to averaged layer 
          !! normal for i_bin.
          write(smcpsi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          psi6_bulk_masked(simbox, particles, (indices == i_bin), layer_normal) 
        end if
      end if

    end do

    !! Append a newline after each row of bins
    write(n_unit, *) ''
    if (write_density) then
      write(density_unit, *) ''
    end if

    if (write_p2) then
      write(p2_unit, *) ''
    end if
    if (write_tau1) then
      write(tau1_unit, *) ''
    end if
    if (write_psi6) then
      write(psi6_unit, *) ''
    end if

    if (write_layerp2) then
      write(layerp2_unit, *) ''
    end if
    if (write_smctau1) then
      write(smctau1_unit, *) ''
    end if
    if (write_smcpsi6) then
      write(smcpsi6_unit, *) ''
    end if
  end do

  call finalize

  contains

  function binvolume(i_bin, maxr, height, binning_direction)
    integer, intent(in) :: i_bin
    real(dp), intent(in) :: maxr
    real(dp), intent(in) :: height
    integer, intent(in) :: binning_direction
    real(dp) :: binvolume
    if (binning_direction == -1) then
      binvolume = 4._dp * atan(1._dp) * height * ((maxr - real(i_bin-1, dp) * &
      bin_width)**2 - max(maxr - real(i_bin, dp) * bin_width, 0._dp)**2)    
    else
      binvolume = 4._dp * atan(1._dp) * height * ((max(real(i_bin, dp) * &
      bin_width, maxr))**2 - (real(i_bin-1, dp) * bin_width)**2)    
    end if    
  end function 

  subroutine finalize
    if(associated(particles)) deallocate(particles)
    if(associated(indices)) deallocate(indices)
    close(configuration_unit)

    if (write_density) then
      close(density_unit)
    end if

    if (write_p2) then
      close(p2_unit)
    end if
    if (write_tau1) then
      close(tau1_unit)
    end if
    if (write_psi6) then
      close(psi6_unit)
    end if

    if (write_layerp2) then
      close(layerp2_unit)
    end if
    if (write_smctau1) then
      close(smctau1_unit)
    end if
    if (write_smcpsi6) then
      close(smcpsi6_unit)
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

