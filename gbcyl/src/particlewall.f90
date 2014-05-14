module particlewall
  use particle, only: particledat, position, orientation 
  use nrtype
  use class_poly_box
  use class_parameterizer
  use class_parameter_writer
  use nrtype
  use ljcylinder
  implicit none
    
  !! @see Micheletti et. al. Journal of Chemical Physics 123, 224705 for 
  !! definitions of these parameters.
  !!
  real(dp), save :: alphaA = 1._dp
  real(dp), save :: alphaB = 1._dp
  real(dp), save :: eps = 1._dp
  real(dp), save :: sig=1._dp
  
  !! The distance of LJ interaction sites from the Gay-Berne particle center
  !! along the unique axis.
  real(dp), save :: LJdist = 1.7_dp 
  logical, save :: isuniformalignment = .false.

  real(dp), save :: epswall_lj = 1._dp
  real(dp), save :: alpha_lj = 1._dp
  real(dp), save :: sigwall_lj = 0.8_dp

  !! General wall properties
  real(dp), save :: wall_density = 1._dp

interface initptwall
  module procedure initptwallparameterizer
end interface

contains 

  subroutine initptwallparameterizer(reader)
    type(parameterizer), intent(in) :: reader
    real(dp) :: Kw_LJ = 5.48819_dp, Kw = 8._dp
    logical :: found
    call getparameter(reader, 'alpha_A', alphaA) 
    call getparameter(reader, 'alpha_B', alphaB) 
    call getparameter(reader, 'LJ_dist', LJdist) 
    call getparameter(reader, 'is_uniform_alignment', isuniformalignment)
    call getparameter(reader, 'sigwall', sig)
    call getparameter(reader, 'alpha_LJ', alpha_lj)
    call getparameter(reader, 'sigwall_LJ', sigwall_lj)
    call getparameter(reader, 'wall_density', wall_density, found)
    if (.not. found) then
      !! Try to read old parameters
      call getparameter(reader, 'Kw', Kw)
      eps = Kw / 8._dp
      call getparameter(reader, 'Kw_LJ', Kw_LJ)
      epswall_lj = Kw_LJ / (Kw * (sigwall_lj/sig)**3) 
    else
      !! Use the new parameters epswall_lj and epswall
      call getparameter(reader, 'epswall_LJ', epswall_lj)
      call getparameter(reader, 'epswall', eps)
    end if 
  end subroutine

  subroutine particlewall_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'General wall parameters')
    call writeparameter(writer, 'wall_density', wall_density)
    !! GB-wall parameters
    call writecomment(writer, 'GB-wall interaction parameters')
    call writeparameter(writer, 'alpha_A', alphaA) 
    call writeparameter(writer, 'alpha_B', alphaB) 
    call writeparameter(writer, 'LJ_dist', LJdist) 
    call writeparameter(writer, 'is_uniform_alignment', isuniformalignment)
    call writeparameter(writer, 'sigwall', sig)
    call writeparameter(writer, 'epswall', eps)
    !! LJ-wall parameters
    call writecomment(writer, 'LJ-wall interaction parameters')
    call writeparameter(writer, 'alpha_LJ', alpha_lj)
    call writeparameter(writer, 'sigwall_LJ', sigwall_lj)
    call writeparameter(writer, 'epswall_LJ', epswall_lj)
  end subroutine


  pure subroutine particlewall_potential(particle, simbox, energy, overlap)
    type(particledat), intent(in) :: particle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    if (particle%rod) then
      call gbwall(particle, simbox, energy, overlap)
    else
      call ljwall(particle, simbox, energy, overlap)
    end if
  end subroutine


  pure function particlewall_force(particle, simbox) result(f)
    type(particledat), intent(in) :: particle
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)
    if (particle%rod) then
      f = gbwall_force(particle, simbox)
    else
      f = ljwall_force(particle, simbox)
    end if
  end function


  !> Calculates the potential energy of a Lennard-Jones (LJ) particle with 
  !! respect to a wall of a cylindrical cavity.
  !! 
  !! @p ljparticle the LJ particle.
  !! @p simbox the simulation box defining the dimensions of the cavity.
  !! @p energy the interaction energy of the wall and the particle.
  !! @p ovrlp tells if the @p ljparticle has penetrated the wall too much.
  !! 
  !! @see D. Micheletti et al. J. Chem. Phys. 123, 224705, 2005.
  !!
  pure subroutine ljwall(ljparticle, simbox, energy, ovrlp)
    implicit none
    type(particledat), intent(in) :: ljparticle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    logical, intent(out) :: ovrlp
    real(dp) :: r
    real(dp) :: Rc
    ovrlp = .false.
    energy = 0._dp
    Rc = getx(simbox)/2._dp 
    r = sqrt(ljparticle%x**2 + ljparticle%y**2)
    if(r >= Rc) then
      ovrlp = .true.
      return
    end if
    energy = ljcylinderpotential(epswall_lj, wall_density, sigwall_lj, alpha_lj, r, Rc)
  end subroutine

  !! Returns the force exerted to a LJ particle by the wall
  pure function ljwall_force(ljparticle, simbox) result(f)
    implicit none
    type(particledat), intent(in) :: ljparticle
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)

    real(dp) :: r
    real(dp) :: Rc
    Rc = getx(simbox)/2._dp 
    r = sqrt(ljparticle%x**2 + ljparticle%y**2)
    f = [ljparticle%x, ljparticle%y, 0._dp]
    f = f / r * ljcylinder_force(epswall_lj, wall_density, sigwall_lj, alpha_lj, r, Rc)
  end function

  !> Calculates the potential energy of a rodlike particle with two embedded 
  !! Lennard-Jones interaction sites with respect to a wall of a cylindrical 
  !! cavity.
  !! 
  !! @p gbparticle the rodlike particle.
  !! @p simbox the simulation box defining the dimensions of the cavity.
  !! @p energy the interaction energy of the wall and the particle.
  !! @p ovrlp tells if the @p gbparticle has penetrated the wall too much.
  !! 
  !! @see D. Micheletti et al. J. Chem. Phys. 123, 224705, 2005.
  !!
  pure subroutine gbwall(gbparticle, simbox, energy,ovrlp)
    implicit none
    type(particledat), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    logical, intent(out) :: ovrlp
    real(dp) :: rsiteA, rsiteB
    real(dp) :: fu, Rc
    ovrlp = .false.
    energy = 0._dp
    Rc = getx(simbox)/2._dp 
    call rArB(gbparticle, rsiteA, rsiteB)
    if(rsiteA >= Rc .or. rsiteB >= Rc) then
      ovrlp = .true.
      return
    end if
    if (isuniformalignment) then
      fu = angular(gbparticle)  
    else 
      fu = 1._dp
    end if
    energy = fu * (ljcylinderpotential(eps, wall_density, sig, alphaA, rsiteA, Rc) + &
    ljcylinderpotential(eps, wall_density, sig, alphaB, rsiteB, Rc))
  end subroutine gbwall

  pure function gbwall_force(gbparticle, simbox) result(f)
    implicit none
    type(particledat), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)
    real(dp) :: fA(3), fB(3)
    real(dp) :: rsiteA, rsiteB
    real(dp) :: fu, Rc
    Rc = getx(simbox)/2._dp 
    call rArB(gbparticle, rsiteA, rsiteB)
    if (isuniformalignment) then
      fu = angular(gbparticle)  
    else 
      fu = 1._dp
    end if

    fA = position(gbparticle) + orientation(gbparticle) * LJdist
    fA(3) = 0._dp
    fA = fA / sqrt(dot_product(fA, fA))

    fB = position(gbparticle) - orientation(gbparticle) * LJdist
    fB(3) = 0._dp
    fB = fB / sqrt(dot_product(fB, fB))

    f = fu * (fA * ljcylinder_force(eps, wall_density, sig, alphaA, rsiteA, Rc) + &
    fB * ljcylinder_force(eps, wall_density, sig, alphaB, rsiteB, Rc))
  end function

  !! Returns the distances of the interaction sites from the axis of the 
  !! cylinder.
  !! 
  !! @p particle the particle of which interaction sites are calcula 
  !! @p ra the distance of site A from the cylinder center
  !! @p rb the distance of site B from the cylinder center
  !!
  pure subroutine rarb(particle, ra, rb)
    implicit none
    intrinsic sqrt
    type(particledat), intent(in) :: particle
    real(dp), intent(out) :: ra, rb
    real(dp) :: xa, ya, xb, yb
    xa=particle%x+LJdist*particle%ux
    ya=particle%y+LJdist*particle%uy
    xb=particle%x-LJdist*particle%ux
    yb=particle%y-LJdist*particle%uy
    ra=sqrt(xa**2 + ya**2)
    rb=sqrt(xb**2 + yb**2)
  end subroutine


  !! Returns the angular dependence of the potential with respect to the 
  !! cylinder axis (z-direction). 
  !!
  !! @p particle the particle to which the potential is calculated. 
  !! 
  pure real(dp) function angular(particle)
    implicit none
    type(particledat), intent(in) :: particle
    angular=(particle%uz)**2
  end function angular

end module particlewall
