program xe_nmr_profile
  use nrtype, only: dp
  use particle, only : particledat
  use orientational_ordering
  use class_poly_box
  use class_factory
  use m_shielding
  use m_quadrupole_coupling
  use m_rank2_tensor, only: rank2_tensor
  use class_parameterizer
  use m_fileunit
  use lj_nmr
  use xfunc_module, only: xfunc => rho
  use histogram
  use utils, only: splitstr, join, fmt_char_dp_array
  implicit none
  type(particledat), allocatable :: particles(:)
  type(particledat), allocatable :: xes(:)
  integer, allocatable :: n_xe(:)
  integer :: io_status
  integer :: i, j, i_bin
  type(poly_box) :: simbox
  type(factory) :: afactory
  type(rank2_tensor) :: gbxe_tensor
  type(rank2_tensor) :: xexe_tensor
  type(rank2_tensor) :: xewall_tensor
  type(rank2_tensor) :: tensor
  real(dp) :: values(3), vectors(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], &
       [3, 3])
  type(parameterizer) :: reader
  integer :: coordinateunit
  logical :: is_wall_on
  character(len=2), parameter :: components(3, 3) = reshape( &
       ["xx", "yx", "zx", "xy", "yy", "zy", "xz", "yx", "zz"], [3, 3])
  character(len=1), parameter :: vi(3) = ["1", "2", "3"]
  procedure(gblj_local), pointer :: gbxe_f
  procedure(ljlj_local), pointer :: xexe_f
  procedure(ljwall_local), pointer :: xewall_f
  real(dp), allocatable :: output(:, :, :)
  integer :: n_tensors 
  real(dp) :: bin_width = -1._dp
  integer :: binning_direction = -1, n_bins = -1
  integer, allocatable :: indices(:)
  logical :: volume_norm = .false.
  logical :: n_xe_norm = .false.
  logical :: is_eigensystem = .false.
  character(len=255) :: configuration_file = '', parameter_file = '', &
       logfile=''
  logical :: q = .false.
  
  !!---- Command line parsing ------------------------------------------------
  character(len=1000) :: cmd_line = ""
  namelist /cmd/ q, configuration_file, parameter_file, bin_width, n_bins, &
       binning_direction, volume_norm, n_xe_norm, is_eigensystem, logfile
  character(len=*), parameter :: options_str = &
       'configuration_file="configurations.txt", ' // &
       'parameter_file="inputparameters.txt", n_bins=100, bin_width=0.1, ' //& 
       '[other options]'

  call get_command(cmd_line)
  call parse_cmdline(cmd_line)

  if (configuration_file == "") then
     write(*, *) 'Give configuration_file!'
     call print_usage
     stop
  else if (parameter_file == "") then
     write(*, *) 'Give parameter_file!'
     call print_usage
     stop
  else if (n_bins < 1) then
     write(*, *) 'Give n_bins > 0!'   
     call print_usage
     stop
  else if (bin_width <= 0._dp) then
     write(*, *) 'Give bin_width > 0!'
     call print_usage
     stop
  end if
  !!--------------------------------------------------------------------------

  if (q) then
    gbxe_f => gbxe_coupling_local
    xexe_f => xexe_coupling_local
    xewall_f => xewall_coupling_local
  else
    gbxe_f => gbxe_shielding_local
    xexe_f => xexe_shielding_local
    xewall_f => xewall_shielding_local
  end if

  if (logfile /= '') then
     reader = new_parameterizer(trim(adjustl(parameter_file)), &
          logfile=logfile)
  else
     reader = new_parameterizer(trim(adjustl(parameter_file)))
  end if

  !! Initialize the modules needed
  call init_shielding(reader)
  call getparameter(reader, "is_wall_on", is_wall_on)
  if (is_wall_on) then
    n_tensors = 4
  else
    n_tensors = 3
  end if
  allocate(output(9, n_tensors, n_bins), n_xe(n_bins))
 
  coordinateunit = fileunit_getfreeunit()
  open(unit=coordinateunit, file=trim(adjustl(configuration_file)), &
    action = 'READ', status = 'OLD')
  if (is_wall_on) then
    write(*, '("#", 1' // fmt_char_dp() //', 44(1X,'// fmt_char_dp() //'))') &
      "x_" // vi, "y_" // vi, "z_" // vi,  "total_" // components, &
      "gbxe_" // components, "xexe_" // components, "xewall_" // components
  else 
    write(*, '("#", 1' // fmt_char_dp() //', 35(1X,'// fmt_char_dp() //'))') &
      "x_" // vi, "y_" // vi, "z_" // vi,  "total_" // components, &
      "gbxe_" // components, "xexe_" // components
  end if
  do  
    call readstate(afactory, coordinateunit, simbox, particles, io_status)  
    if (io_status < 0) exit
    if (is_eigensystem) then
       call eigens(pack(particles, particles%rod), count(particles%rod), &
            values, vectors)
       call cycle_largest_to_3rd(values, vectors)
    end if
    !! Allocate memory for indices
    if (.not. allocated(indices)) then
      allocate(indices(size(particles)))
    else if (size(indices) /= size(particles)) then
      if (allocated(indices)) deallocate(indices)
      allocate(indices(size(particles)))
    end if  
    !! Solve indices

    if (binning_direction == -1) then
      call bin_indices(particles, size(particles), xfunc, bin_width, indices, &
      getx(simbox) / 2._dp, binning_direction)
    else 
      call bin_indices(particles, size(particles), xfunc, bin_width, indices)
    end if

    n_xe = 0
    do i_bin = 1, n_bins
      n_xe(i_bin) = count((.not. particles%rod) .and. indices == i_bin)
            
      if (.not. allocated(xes)) then
        allocate(xes(n_xe(i_bin)))
      else if (size(xes) /= n_xe(i_bin)) then
        if (allocated(xes)) deallocate(xes)
        allocate(xes(n_xe(i_bin)))
      end if
      xes = pack(particles, (.not. particles%rod) .and. indices == i_bin)

      !! Accumulate Xe-GB shielding
      gbxe_tensor = new_tensor(axes = vectors)
      do i = 1, n_xe(i_bin)
         do j = 1, size(particles)
            if (particles(j)%rod) gbxe_tensor = gbxe_tensor + &
                 rotate(gblj_tensor(gbxe_f, simbox, particles(j), xes(i)), &
                 gbxe_tensor)
         end do
      end do

      !! Accumulate Xe-Xe shielding
      xexe_tensor = new_tensor(axes = vectors)
      do i = 1, n_xe(i_bin) - 1
         !! Iterate through Xenons in the same bin
         do j = i + 1, n_xe(i_bin)
            xexe_tensor = xexe_tensor + rotate(ljlj_tensor(&
                 xexe_f, simbox, xes(i), xes(j)), xexe_tensor)       
         end do
      end do

      do i = 1, n_xe(i_bin)
         !! Iterate through Xenons in other bins:
         do j = 1, size(particles)
            if ((.not. (indices(j) == i_bin)) .and. (.not. particles(j)%rod)) &
               then
               !! Should we only add half of the contribution to this bin to 
               !! avoid double-counting? No, because we are not calculating 
               !! the  total "NMR-energy" but an average of the interactions 
               !! "felt" by single atoms
               xexe_tensor = xexe_tensor + rotate(ljlj_tensor(&
                    xexe_f, simbox, xes(i), particles(j)), xexe_tensor)       
            end if
         end do
      end do

      tensor = xexe_tensor + gbxe_tensor

      if (is_wall_on) then
         !! Compute wall contribution 
         xewall_tensor = new_tensor(axes = vectors)
         do i = 1, n_xe(i_bin)
            xewall_tensor = xewall_tensor + rotate(ljwall_tensor( &
                 xewall_f, simbox, xes(i)), xewall_tensor)
         end do
         tensor = tensor + xewall_tensor
         output(:, 4, i_bin) = reshape(xewall_tensor%components, &
              shape(output(:, 4, i_bin)))
      end if
      output(:, 1, i_bin) = reshape(tensor%components, &
           shape(output(:, 1, i_bin)))
      output(:, 2, i_bin) = reshape(gbxe_tensor%components, &
           shape(output(:, 2, i_bin)))
      output(:, 3, i_bin) = reshape(xexe_tensor%components, &
           shape(output(:, 3, i_bin)))

      !! If averaging over xenons in each configurations is desired:
      if (n_xe_norm) output(:, :, i_bin) = output(:, :, i_bin) / &
           real(n_xe(i_bin), dp)

      !! Normalize by bin volume. So effectively we are calculating a 
      !! shielding/coupling density.
      if (volume_norm) output(:, :, i_bin) = output(:, :, i_bin) / &
           binvolume(i_bin, getx(simbox)/2._dp, getz(simbox), &
           binning_direction) !* volume(simbox) / count(.not. particles%rod)

    end do !! end going through bins for this snapshot
    write(*, fmt = fmt_char_dp_array(size(output) + size(vectors) + n_bins)) &
         vectors, (n_xe(i_bin), output(:, :, i_bin), i_bin = 1, n_bins)
  end do

  call fileunit_closeallopen()

contains 

include 'parse_cmdline.inc'

subroutine print_help
  call print_usage
  write(*, *) 'This program computes the shielding/quadrupole coupling tensor'
  write(*, *) 'for a 129/131Xe atom dissolved in a Gay-Berne(4.4, 20, 1, 1)'
  write(*, *) 'liquid crystal confined to a cylindrical cavity. The Xe'
  write(*, *) 'shielding/quadrupole coupling is calculated by dividing the'
  write(*, *) 'system into cylindrical shells (bins) and accumulating (or '
  write(*, *) 'averaging if normalization is used -- see options below -- )'
  write(*, *) 'the tensors in each bin.'
  write(*, *) ''
  write(*, *) "The options for the program are listed below." 
  write(*, *) ""
  write(*, *) "configuration_file=filename"
  write(*, *) "    molecular coordinates are given in this file"
  write(*, *) ""
  write(*, *) "parameter_file=filename"
  write(*, *) "    input parameters of the simulation are given in this file"
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
  write(*, *) 'q=F' 
  write(*, *) '    when q=T, quadrupole coupling is computed in place of'
  write(*, *) '    shielding.'
  write(*, *) ''
  write(*, *) 'volume_norm=F'
  write(*, *) '    when volume_norm=T, shielding is normalized with bin volume'
  write(*, *) ''
  write(*, *) 'n_xe_norm=F'
  write(*, *) '    when n_xe_norm, shielding is normalized with the number of'
  write(*, *) '    xe in each bin. Note that this can result in NaN in the'
  write(*, *) '    output'
  write(*, *) ''
  write(*, *) 'is_eigensystem=F'
  write(*, *) '    when is_eigensystem=T, the shielding/quadrupole coupling is'
  write(*, *) '    computed in the principal axis system of the orientational'
  write(*, *) '    ordering tensor, so that the eigenvector corresponding to'
  write(*, *) '    the largest eigenvalue of the orientational ordering tensor'
  write(*, *) '    is the 3rd axis.'
  write(*, *) ''
  write(*, *) 'logfile=filename'
  write(*, *) '    log about parameter reading. Useful for reducing clutter in'
  write(*, *) '    the output of the program.'
  write(*, *) ''
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

end program


