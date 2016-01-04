program xe_nmr
  use num_kind, only: dp
  use m_particledat, only : particledat
  use orientational_ordering, only: orientation_parameter, eigens, &
       cycle_largest_to_3rd
  use class_poly_box, only: poly_box
  use m_particle_factory, only: factory, factory_readstate
  use m_shielding, only: gbxe_shielding_local, xexe_shielding_local, &
       xewall_shielding_local, init_shielding
  use m_quadrupole_coupling, only: gbxe_coupling_local, xexe_coupling_local, &
       xewall_coupling_local
  use m_rank2_tensor
  use class_parameterizer
  use m_fileunit
  use lj_nmr, only: gblj_tensor, ljlj_tensor, ljwall_tensor, gblj_local, &
       ljlj_local, ljwall_local
  use utils, only: splitstr, join, fmt_char_dp
  implicit none
  type(particledat), allocatable :: particles(:)
  type(particledat), allocatable :: xes(:)
  integer :: n_xe
  integer :: io_status
  integer :: i, j
  type(poly_box) :: simbox
  type(factory) :: afactory
  type(rank2_tensor) :: gbxe_tensor
  type(rank2_tensor) :: xexe_tensor
  type(rank2_tensor) :: xewall_tensor
  type(rank2_tensor) :: tensor
  real(dp) :: values(3), vectors(3, 3)
  type(parameterizer) :: reader
  character(4) :: idchar = ""
  integer :: coordinateunit
  logical :: is_wall_on
  character(len=2), parameter :: indices(3, 3) = reshape(["xx", "yx", "zx", &
    "xy", "yy", "zy", "xz", "yx", "zz"], [3, 3])
  character(len=1), parameter :: vi(3) = ["1", "2", "3"]
  procedure(gblj_local), pointer :: gbxe_f
  procedure(ljlj_local), pointer :: xexe_f
  procedure(ljwall_local), pointer :: xewall_f

  !!---- Command line parsing ------------------------------------------------
  character(len=255) :: cmd_line = ""
  character(len=80), allocatable :: args(:)
  character(len=255) :: message
  logical :: q = .false.
  namelist /cmd/ q, idchar
  integer :: ios
  call get_command(cmd_line)
  !! Remove program name from cmd_line
  call splitstr(cmd_line, " ", args)
  call join(args(2:), "", cmd_line)
  cmd_line = "&cmd " // trim(adjustl(cmd_line)) // " /" 
  read(cmd_line, NML = cmd, iostat = ios, iomsg = message)
  if (ios /= 0) then
    write(*, *) "ERROR: ", ios
    write(*, *) message
  end if
  if (idchar == "") stop 'Usage: ./xe_nmr [q=T] idchar="0"'
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

  reader = new_parameterizer('inputparameters.'//trim(adjustl(idchar)), &
    logfile = "xe_shielding_log."//trim(adjustl(idchar)))
  !! Initialize the modules needed
  call init_shielding(reader)
  call getparameter(reader, "is_wall_on", is_wall_on)
  coordinateunit = fileunit_getfreeunit()
  open(unit=coordinateunit, file='configurations.' //trim(adjustl(idchar)), &
    action = 'READ', status = 'OLD')
  if (is_wall_on) then
    write(*, '("#", 1' // fmt_char_dp() //', 44(1X,'// fmt_char_dp() //'))') &
      "x_" // vi, "y_" // vi, "z_" // vi,  "total_" // indices, &
      "gbxe_" // indices, "xexe_" // indices, "xewall_" // indices
  else 
    write(*, '("#", 1' // fmt_char_dp() //', 35(1X,'// fmt_char_dp() //'))') &
      "x_" // vi, "y_" // vi, "z_" // vi,  "total_" // indices, &
      "gbxe_" // indices, "xexe_" // indices
  end if
  do  
    call factory_readstate(afactory, coordinateunit, simbox, particles, io_status)  
    if (io_status < 0) exit
    call eigens(pack(particles, particles%rod), count(particles%rod), &
      values, vectors)
    call cycle_largest_to_3rd(values, vectors)

    n_xe = count(.not. particles%rod)
    if (.not. allocated(xes) .or. n_xe /= size(xes)) allocate(xes(n_xe))
    xes = pack(particles, .not. particles%rod)


    !! Average Xe-GB shielding
    gbxe_tensor = new_tensor(axes = vectors)
    do i = 1, n_xe
      do j = 1, size(particles)
        if (particles(j)%rod) gbxe_tensor = gbxe_tensor + rotate(gblj_tensor(&
          gbxe_f, simbox, particles(j), xes(i)), gbxe_tensor)
      end do
    end do
    gbxe_tensor = gbxe_tensor / n_xe

    !! Average Xe-Xe shielding
    xexe_tensor = new_tensor(axes = vectors)
    do i = 1, n_xe - 1
      do j = i + 1, n_xe
        xexe_tensor = xexe_tensor + rotate(ljlj_tensor(&
          xexe_f, simbox, xes(i), xes(j)), xexe_tensor)       
      end do
    end do
    xexe_tensor = xexe_tensor / n_xe

    tensor = xexe_tensor + gbxe_tensor

    if (is_wall_on) then
      !! Compute wall contribution 
      xewall_tensor = new_tensor(axes = vectors)
      do i = 1, n_xe
        xewall_tensor = xewall_tensor + rotate(ljwall_tensor( &
          xewall_f, simbox, xes(i)), xewall_tensor)
      end do
      xewall_tensor = xewall_tensor / n_xe
      tensor = tensor + xewall_tensor
      !! Print instantaneous axes (X, Y, Z), the total shielding tensor and 
      !! individual contributions from GB and Xe particles and the wall 
      !! (if present)
      write(*, '( 45('// fmt_char_dp() //',1X))') vectors(1:3, 1), &
        vectors(1:3, 2), vectors(1:3, 3), tensor%components, &
        gbxe_tensor%components, xexe_tensor%components, &
        xewall_tensor%components
    else
      write(*, '( 36('// fmt_char_dp() //',1X))') vectors(1:3, 1), &
        vectors(1:3, 2), vectors(1:3, 3), tensor%components, &
        gbxe_tensor%components, xexe_tensor%components
    end if      
  end do
end program


