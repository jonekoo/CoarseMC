!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Divides the system in to cylindrical bins and calculates for each bin
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
!! These latter ones have the average layer normal as their reference 
!! direction.
!! 
!! For options and their descriptions run with 
!! ./cylhist --help 
!! 
program cylhist
  use nrtype
  use particle
  use orientational_ordering
  use tau1_module
  use tau1_negative
  use histogram
  use utils
  use class_factory
  use class_poly_box
  use psi6_module
  use layernormal
  use m_fileunit
  implicit none

  type(particledat), allocatable :: particles(:), temp(:)
  integer :: io_status
  integer :: n_particles
  
  real(dp) :: bin_width = 0.5_dp
  integer :: binning_direction = 1

  real(dp) :: bin_vol

  !! If bins with uniform volume are desired
  logical :: uniform_volume = .false.
  !real(dp), allocatable :: boundaries(:)
  !real(dp), allocatable :: weighted_midpoints(:)

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
  character(len=200) :: tau1sqr_file = ''
  character(len=200) :: r_file = ''

  integer :: configuration_unit
  integer :: psi6_unit 
  integer :: layerp2_unit
  integer :: smctau1_unit
  integer :: smcpsi6_unit
  integer :: p2_unit
  integer :: tau1_unit
  integer :: density_unit
  integer :: n_unit
  integer :: tau1sqr_unit
  integer :: r_unit

  character(80) :: particletype = 'gb'

  namelist /cmd/ bin_width, cutoff, psi6_file, smcpsi6_file, &
       smctau1_file, layerp2_file, binning_direction, p2_file, tau1_file, &
       density_file, n_file, tau1sqr_file, configuration_file, n_bins, particletype, &
       uniform_volume, r_file
       
  integer, dimension(:), pointer :: indices 
  integer :: n_bins = 100
  integer :: i_bin
  integer :: allocstat
  integer :: n_bin_particles
  real(dp) :: tensor(3,3)
  real(dp) :: layer_p2
  real(dp) :: layer_normal(3)
  real(sp) :: bin_layer_normal(3) 
  real(sp) :: tau1_value
  real(sp) :: smctau1_value
  real(sp) :: layer_distance
  real(sp) :: global_layer_distance
  complex(spc) :: global_tau1c, tau1c
  real(sp) :: smclayer_distance
  real(sp), dimension(3), parameter :: z_axis=(/0._sp, 0._sp, 1._sp/)
  !!real(sp) :: tau1sqr_sum = 0.0

  !! Parameters to control which observables are calculated and written.
  logical :: write_psi6 = .false.
  logical :: write_smcpsi6 = .false.
  logical :: write_smctau1 = .false.
  logical :: write_layerp2 = .false.
  logical :: write_p2 = .false.
  logical :: write_density = .false.
  logical :: write_tau1 = .false.
  logical :: write_particle_count = .false.
  logical :: write_tau1sqr = .false.
  logical :: write_r = .false.

  type(poly_box) :: simbox
  type(factory) :: afactory
  integer :: i_conf

  character(len=*), parameter :: options_str = &
       'configuration_file=configurations.0, n_bins=100, bin_width=0.1, ' // &
       'n_file, particletype="gb", [other_options]'
  character(len=1000) :: cmd_line = ""

  call get_command(cmd_line)
  call parse_cmdline(cmd_line)

  write(*, *) 'cylhist2 called with parameters:'
  write(*, nml=cmd)

  if (configuration_file == "") then
     write(*, *) 'Give configuration_file!'
     call print_usage()
     stop
  else if (n_bins < 1) then
     write(*, *) 'Give n_bins >= 1 !'
     call print_usage()
     stop
  else if(bin_width <= 0._dp) then
     write(*, *) 'Give bin_width > 0 !'
     call print_usage() 
     stop
  else if (trim(n_file) == "") then
     write(*, *) 'Give n_file!'
     call print_usage
     stop 
  end if

  call open_write_unit(n_file, n_unit)
  write_particle_count = .true.

  !! Open configuration_file
  if (trim(configuration_file) /= '') then
    configuration_unit=fileunit_getfreeunit()
    open(unit=configuration_unit, file=configuration_file, action='READ', &
         status='OLD')
  else 
    write(*, *) 'Set configuration_file!' 
    stop
  end if


  if (trim(r_file) /= '') then
     call open_write_unit(r_file, r_unit)
     write_r = .true.
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
  if (trim(tau1sqr_file) /= '') then
    call open_write_unit(tau1sqr_file, tau1sqr_unit)
    write_tau1sqr = .true.
  end if

  call init_histogram(bin_width, n_bins, uniform_volume)
  i_conf = 0
  do  
    call factory_readstate(afactory, configuration_unit, simbox, temp, io_status)
    if (io_status /= 0) then
      exit
    end if

    if (allocated(particles)) deallocate(particles)
    if (particletype == 'lj') then
      n_particles = count(.not. temp%rod)
      allocate(particles(n_particles))
      particles(1:n_particles) = pack(temp, .not. temp%rod)
      if(allocated(temp)) deallocate(temp)
    else if (trim(adjustl(particletype)) == 'gb') then
      n_particles = count(temp%rod)
      allocate(particles(n_particles))
      particles(1:n_particles) = pack(temp, temp%rod)
      if(allocated(temp)) deallocate(temp)
    else 
      n_particles = size(temp)
      call move_alloc(temp, particles)
    end if
    i_conf = i_conf + 1
    if (i_conf > 1) then
       !! Append a newline after each row of bins
       write(n_unit, '(/)', ADVANCE='NO')
       if (write_density) then
          write(density_unit, '(/)', ADVANCE='NO')
       end if

       if (write_p2) then
          write(p2_unit, '(/)', ADVANCE='NO')
       end if
       if (write_tau1) then
          write(tau1_unit, '(/)', ADVANCE='NO')
       end if
       if (write_tau1sqr) then
          write(tau1sqr_unit, '(/)', ADVANCE='NO')
       end if
       if (write_psi6) then
          write(psi6_unit, '(/)', ADVANCE='NO')
       end if
       
       if (write_layerp2) then
          write(layerp2_unit, '(/)', ADVANCE='NO')
       end if
       if (write_smctau1) then
          write(smctau1_unit, '(/)', ADVANCE='NO')
       end if
       if (write_smcpsi6) then
          write(smcpsi6_unit, '(/)', ADVANCE='NO')
       end if
       if (write_r) then
          write(r_unit, '(/)', ADVANCE='NO')
       end if
    end if

    if (write_r) then
       write(r_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') avg_r(particles)
    end if

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
 
    if (write_tau1sqr) then
      call tau1_routine(particles, z_axis, tau1_value, global_layer_distance)
      global_tau1c = tau1_c(real(particles%z, sp), global_layer_distance)
      write(tau1sqr_unit, '(4('//fmt_char_dp()//',1X))', ADVANCE='NO') &
           real(conjg(global_tau1c) * global_tau1c), size(particles), &
           real(global_tau1c), aimag(global_tau1c)
      !! The written value should be equal to the global tau1 with respect to
      !! z-axis.
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
      call tau1_routine(particles, real(layer_normal, sp), smctau1_value, &
           smclayer_distance)
      write(smctau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
           smctau1_value
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

    call bin_indices(particles, indices)

    do i_bin = 1, n_bins
      n_bin_particles = count(indices == i_bin)

      write(n_unit, '('//fmt_char_int()//',1X)', ADVANCE='NO') n_bin_particles

      if (write_r) then
         write(r_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
              avg_r(pack(particles, indices == i_bin))
      end if

      if (write_density) then
         !! Calculate and write density profile
         bin_vol = getz(simbox) * bin_base_area(i_bin)
         write(density_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
              real(n_bin_particles, dp) / (bin_base_area(i_bin) * getz(simbox))
      end if

      if (write_p2) then
        if(n_bin_particles == 0) then    
          write(p2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 'NaN'
        else
          !! Calculate and write orientation parameter profile
          write(p2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
               orientation_parameter_v(pack(particles(1:n_particles), &
               indices==i_bin), real(z_axis, dp))
        end if
      end if

      if (write_tau1) then
        if(n_bin_particles == 0) then    
          write(tau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 'NaN'
        else
          !! Calculate and write tau1 profile        
          write(tau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          tau1(pack(particles(1:n_particles), indices==i_bin), z_axis, &
          layer_distance)
        end if
      end if

      if (write_tau1sqr) then
        if(n_bin_particles == 0) then    
          write(tau1sqr_unit, '(4('//fmt_char_dp()//',1X))', ADVANCE='NO')&
               'NaN', 'NaN', 'NaN', 'NaN'
        else
          !! Calculate and write tau1 profile        
          !! Should it be sqrt(real()) or abs()? When the complex numbers are
          !! equal it does not matter, but otherwise it makes a difference.
          tau1c = tau1_c(pack(real(particles%z, sp), indices==i_bin), &
               global_layer_distance)
          write(tau1sqr_unit, '(4('//fmt_char_dp()//',1X))', ADVANCE='NO') &
               real(conjg(global_tau1c) * tau1c), n_bin_particles, &
               real(tau1c), aimag(tau1c)
        end if
      end if

      if (write_psi6) then
        if(n_bin_particles == 0) then         
          write(psi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 'NaN'
        else
          !! Calculate and write psi6 profile
          write(psi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          psi6_bulk_masked(simbox, particles, (indices == i_bin), real(z_axis, dp)) 
        end if
      end if

      if (write_layerp2) then
        if(n_bin_particles == 0) then    
          write(layerp2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 'NaN'
        else
          !! Write the orientation parameter for layer normals.
          write(layerp2_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          dot_product(bin_layer_normal, layer_normal)
        end if
      end if

      if (write_smctau1) then
        if(n_bin_particles == 0) then    
          write(smctau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 'NaN'
        else
          !! Calculate and write smctau1 profile        
          write(smctau1_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          tau1(pack(particles(1:n_particles), indices==i_bin), &
          real(layer_normal, sp), smclayer_distance)
        end if
      end if
    
      if (write_smcpsi6) then
        if(n_bin_particles == 0) then    
          write(smcpsi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') 'NaN'
        else
          !! Calculate and write psi6 profile with respect to averaged layer 
          !! normal for i_bin.
          write(smcpsi6_unit, '('//fmt_char_dp()//',1X)', ADVANCE='NO') &
          psi6_bulk_masked(simbox, particles, (indices == i_bin), layer_normal)
        end if
      end if

    end do
    !! debug stuff
    !!write(*, *) 'tau1sqr_sum = ', tau1sqr_sum,  'global_tau1sqr = ', &
    !!     conjg(global_tau1c) * global_tau1c, 'tau1_sqr = ', tau1_value**2 
    !!tau1sqr_sum = 0.0
  end do

  call finalize

contains

include 'parse_cmdline.inc'

subroutine print_help
  call print_usage
  write(*, *) "This is a program to compute structural and order parameters"
  write(*, *) "for a Gay-Berne liquid crystal (and LJ particles) confined to"
  write(*, *) "a cylindrical cavity. The properties are calculated by dividing"
  write(*, *) "the molecules to cylindrical shells (bins) of thickness"
  write(*, *) "bin_width. This gives us information about them as a function of"
  write(*, *) "distance to the wall or to the cavity center. Results are"
  write(*, *) "printed to separate outputfiles for each quantity. In the"
  write(*, *) "output files one line corresponds to one snapshot of molecules"
  write(*, *) "one record is an average over the molecules in one bin. Order"
  write(*, *) "is determined with option binning_direction."
  write(*, *) ""
  write(*, *) ""
  write(*, *) "The options for the program are listed below." 
  write(*, *) ""
  write(*, *) "configuration_file"
  write(*, *) "    molecular coordinates are given in this file"
  write(*, *) ""
  write(*, *) "n_file=filename"
  write(*, *) "    file to write number of molecules in each bin."
  write(*, *) ""
  write(*, *) "bin_width=number"
  write(*, *) "    thickness of the cylindrical shell"
  write(*, *) ""
  write(*, *) "binning_direction=1 or -1"
  write(*, *) "    if -1 the first bin is the outermost and if 1 the first"
  write(*, *) "    bin is the innermost."
  write(*, *) ""
  write(*, *) "n_bins=integer"
  write(*, *) "    gives the number of cylindrical shells into which the"
  write(*, *) "    system is divided. Restrictions: n_bins > 0 and" 
  write(*, *) "    n_bins * bin_width <= system radius."
  write(*, *) ""
  write(*, *) "particletype='gb'"
  write(*, *) "    the particle type which are subjected to computations. If"
  write(*, *) "    one wants to compute the density_profile for LJ particles"
  write(*, *) "    he/she should select 'lj' default value is 'gb'"
  write(*, *) ""
  write(*, *) "All of the following are options which set an outputfile."
  write(*, *) "the name of the option is related to the calculated quantity."
  write(*, *) "None of them are necessary to run the program and should be"
  write(*, *) "given only if the calculation of the respected quantity is"
  write(*, *) "desired. All filenames must be different!"
  write(*, *) ""
  write(*, *) "p2_file=filename"
  write(*, *) "    orientational ordering parameter"
  write(*, *) ""
  write(*, *) "tau1_file=filename"
  write(*, *) "    1D translational order parameter"
  write(*, *) ""
  write(*, *) "psi6_file=filename"
  write(*, *) "    bond-orientational order parameter"
  write(*, *) ""
  write(*, *) "density_file=filename"
  write(*, *) "    number density"
  write(*, *) ""
  write(*, *) "layerp2_file=filename"
  write(*, *) "    order parameter for molecular layers."
  write(*, *) ""
  write(*, *) "smctau1_file=filename"
  write(*, *) "    translational order for possibly tilted smectic-C like" 
  write(*, *) "    layers."
  write(*, *) ""
  write(*, *) "smcpsi6_file=filename"
  write(*, *) "    bond-orientational order in tilted layers"
  write(*, *) ""
  write(*, *) "Other options:"
  write(*, *) ""
  write(*, *) "cutoff=number"
  write(*, *) "    determines a cutoff radius to the molecules which are"
  write(*, *) "    included in the local layer for the layerp2."
  write(*, *) ""
end subroutine


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
end function binvolume


function avg_r(particles)
  type(particledat), intent(in) :: particles(:)
  real(dp) :: avg_r
  integer :: i
  avg_r = 0
  do i = 1, size(particles)
     avg_r = avg_r + sqrt(particles(i)%x**2 + particles(i)%y**2)
  end do
  avg_r = avg_r / size(particles)
end function


subroutine finalize
  if(associated(indices)) deallocate(indices)
  close(configuration_unit)
  
  if (write_r) close(r_unit)

  if (write_density) then
     close(density_unit)
  end if
  
  if (write_p2) then
     close(p2_unit)
  end if
  if (write_tau1) then
     close(tau1_unit)
  end if
  if (write_tau1sqr) then
     close(tau1sqr_unit)
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
end subroutine finalize


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


end program cylhist

