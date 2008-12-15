program gbcylinder
  use nrtype
  use mcstep, only : step, updatemaxvalues, initmcstep, &
                     & mcstep_write_state => write_module_state
  use energy, only : totenergy, totwallprtclV, totpairV
  use cylinder, only : initcylinder, getradius, getHeight
  use verlet, only : initvlist, freevlist
  use particlewall, only : initptwall, rArB
  use gbgb, only : initgbgb
  use particle, only : initParticle, &
                     & particle_write_state => write_module_state, &
                     & particledat
  use io, only : readState, io_writeState => writestate, povout, ReadParams, &
               & initOutput 
  use mtmod, only: sgrnd
  implicit none

  integer :: anchor 
  integer, parameter :: eunit = 10
  integer :: opening_status
  character(LEN = *), parameter :: efile='energy.txt'
  character(len = 50) :: statefile 
  type(particledat), dimension(:), pointer :: particleArray
  integer :: N, Nrelax, Nprod, Nratio, seed, i, vtype
  real(dp) :: radius, height, Kw
  real(dp) :: T, pres, epses, eps0, rsphere, spmyy, epsphere, sigma0, siges
  real(dp) :: totE = 0.0
  real(dp) ::  maxangle = 0.170, maxtrans = 0.156
  logical, target :: ol = .false.  
  logical, pointer :: overlap
  integer :: debug, allign

  overlap => ol 
  call ReadParams(statefile, Nrelax, Nprod, Nratio, T, pres, anchor, vtype, &
                & Kw, seed, epses, eps0, rsphere, spmyy, epsphere, sigma0, & 
                & siges, allign, debug)  
  call readstate(statefile, particleArray, radius, height)
  !! Initialize modules. 
  call initptwall(anchor, Kw)
  call sgrnd(seed)  
  call initgbgb()
  call initOutput()
  call initcylinder(radius, height, T, pres)
  call initmcstep(vtype)
  call initvlist(particleArray)
  call initParticle(maxtrans, maxangle)

  open(eunit, FILE = efile, status = 'replace', position = 'append', &
     & iostat = opening_status)
  if ( opening_status /= 0) then
    write (*,*) 'Trying to open ', efile,' resulted in an error.'
    write (*,*) 'Program will end.'
    stop
  end if
  write(eunit, *) 'MC-sweep  Total energy'

  do i = 1, (Nrelax + Nprod) 
    if (mod(i, Nratio)==0) then
      radius = getRadius()
      height = getHeight()
      call totEnergy(particleArray, overlap, totE)
      write (eunit,*) i, totE
      call io_writestate(radius, height, particleArray)
      if (i .le. Nrelax) then !! Equilibration only
         N = size(particleArray) 
         call updateMaxValues(N, Nratio)
      end if
    end if
    call step(particleArray)
  end do
   
  !! Kirjoitetaan hiukkasten loppupaikat muistiin. 
  call povout(particleArray, radius, height)
  close(eunit)
  call freevlist()
  if (associated(particleArray)) deallocate(particleArray)
  write (*, *) 'End of program.'
  
end program gbcylinder
