!! Usage: program < configuration_file
!!
!! Program reads the given file configuration by configuration. a For each 
!! configuration it calls the routines calculate and print which redirects
!! calls to analysis routines and prints the results to files. 
!!
program analysis
  use state_reader, only: read_configuration
  use nrtype, only: dp
  use io, only: ReadParams
  use particlewall, only: initptwall
  use gayberne, only: gayberne_init => init
  use particle, only: particledat
  use particle_mover, only: initParticle
  use cylinder, only: initcylinder
  use verlet, only: initvlist, freevlist
  use energy, only: total_energy
  implicit none
  type(particledat), dimension(:), pointer :: particles
  integer :: io_status
  integer :: n_particles
  !! Initialize energy calculation.  
    character(len = 78) :: statefile 
    real(dp) :: temperature
    real(dp) :: pressure
    integer :: anchor
    integer :: volume_scaling_type
    real(dp) :: Kw
    integer :: seed
    real(dp) :: kappa_epsilon_sgb
    real(dp) :: epsilon_0_sgb
    real(dp) :: r_sphere
    real(dp) :: spmyy
    real(dp) :: epsilon_ss
    real(dp) :: sigma_0_sgb
    real(dp) :: kappa_sigma_sgb
    integer :: allign
    integer :: debug 
    real(dp) :: radius
    real(dp) :: height
    integer :: adjusttype
    real(dp) :: move_ratio
    real(dp) :: scaling_ratio
    real(dp) :: max_translation
    real(dp) :: max_rotation
    real(dp), parameter :: kappa_sigma = 4.4
    real(dp), parameter :: kappa_epsilon = 20.0
    real(dp), parameter :: mu = 1.0
    real(dp), parameter :: nu = 1.0
    real(dp), parameter :: sigma_0 = 1.0
    real(dp), parameter :: epsilon_0 = 1.0
    real(dp), dimension(3), parameter :: magnetic_field_direction = &
      & (/0.0, 0.0, 1.0/)
    real(dp), parameter :: magnetic_field_teslas = 11.74
    logical, parameter :: is_magnet_on = .false.
    integer :: n_equilibration_sweeps_
    integer :: n_production_sweeps_
    integer :: production_period_
    logical :: overlap
    real(dp) :: total_e
    call ReadParams(statefile, n_equilibration_sweeps_, n_production_sweeps_, &
      & production_period_, temperature, pressure, anchor, & 
      & volume_scaling_type, Kw, seed, kappa_epsilon_sgb, epsilon_0_sgb, & 
      & r_sphere, spmyy, epsilon_ss, sigma_0_sgb, kappa_sigma_sgb, allign, & 
      & debug, adjusttype, move_ratio, scaling_ratio, max_translation, &
      & max_rotation)
    call initptwall(anchor, Kw)
    call gayberne_init(kappa_sigma, kappa_epsilon, mu, nu, sigma_0, epsilon_0)
  do  
    call read_configuration(particles, n_particles, radius, height, io_status)
    if (io_status < 0) then
      exit
    else
      call initcylinder(radius, height)
      call initvlist(particles, n_particles)
      call total_energy(particles, n_particles, overlap, total_e)
      write(*,*) total_e 
      call freevlist
    end if
  end do
  if(associated(particles)) deallocate(particles)
end program analysis


