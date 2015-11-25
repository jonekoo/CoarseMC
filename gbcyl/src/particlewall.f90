!> Module for computing particle-wall interactions for a Lennard-Jones
!! (LJ)  or Gay-Berne (GB) particle in a cylindrical cavity. The cavity
!! walls are smooth and consist of evenly distributed LJ particles. A
!! GB particle interacts with the wall through two embedded LJ sites in
!! the particle.
!! 
!! @see Micheletti et. al. Journal of Chemical Physics 123, 224705.
!!
module particlewall
  use particle, only: particledat, position, orientation 
  use num_kind
  use class_poly_box
  use class_parameterizer
  use class_parameter_writer
  use ljcylinder
  use json_module
  implicit none
    
  !> The relative strength of attractive term as compared to the
  !! repulsive term for the LJ site A of the GB particle. 
  real(dp), save :: alpha_a = 1.

  !> The relative strength of attractive term as compared to the
  !! repulsive term for the LJ site B of the GB particle. 
  real(dp), save :: alpha_b = 1.

  !> The well-depth parameter of the interaction between the wall and the
  !! LJ site.
  real(dp), save :: eps = 1.

  !> The range parameter of the interaction between the LJ site and the
  !! wall.
  real(dp), save :: sig = 1.
  
  !> The distance of LJ interaction sites from the Gay-Berne particle
  !! center along the unique axis.
  real(dp), save :: LJdist = 1.7 

  !> If true, the wall interaction favors uniform alignment of GB 
  !! particles along the z-axis of the simulation box.
  logical, save :: isuniformalignment = .false.

  !> Well-depth parameter for the interaction of a LJ particle with the
  !! wall.
  real(dp), save :: epswall_lj = 1.

  !> The strength of the attractive term with respect to the repulsive
  !! term for the interaction of a LJ particle with the wall.
  real(dp), save :: alpha_lj = 1.
  
  !> Range parameter for the interaction of a LJ particle with the
  !! wall.
  real(dp), save :: sigwall_lj = 0.8

  !> Number density of virtual LJ particles in the wall.
  real(dp), save :: wall_density = 1.

  !> Initializes the module.
  interface particlewall_init
     module procedure particlewall_init1, particlewall_from_json
  end interface

contains 

  !> Initializes the module with parameters read by @p reader.
  subroutine particlewall_init1(reader)
    type(parameterizer), intent(in) :: reader
    real(dp) :: Kw_LJ = 5.48819_dp, Kw = 8._dp
    logical :: found
    call getparameter(reader, 'alpha_A', alpha_a) 
    call getparameter(reader, 'alpha_B', alpha_b) 
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

  !> Initializes the module with parameters read by @p reader.
  subroutine particlewall_from_json(json_val)
    type(json_value), pointer, intent(in) :: json_val
    real(dp) :: Kw_LJ = 5.48819_dp, Kw = 8._dp
    logical :: found
    call json_get(json_val, 'alpha_A', alpha_a) 
    call json_get(json_val, 'alpha_B', alpha_b) 
    call json_get(json_val, 'LJ_dist', LJdist) 
    call json_get(json_val, 'is_uniform_alignment', isuniformalignment)
    call json_get(json_val, 'sigwall', sig)
    call json_get(json_val, 'alpha_LJ', alpha_lj)
    call json_get(json_val, 'sigwall_LJ', sigwall_lj)
    call json_get(json_val, 'wall_density', wall_density, found)
    if (.not. found) then
      !! Try to read old parameters
      call json_get(json_val, 'Kw', Kw)
      eps = Kw / 8._dp
      call json_get(json_val, 'Kw_LJ', Kw_LJ)
      epswall_lj = Kw_LJ / (Kw * (sigwall_lj/sig)**3) 
    else
      !! Use the new parameters epswall_lj and epswall
      call json_get(json_val, 'epswall_LJ', epswall_lj)
      call json_get(json_val, 'epswall', eps)
    end if 
  end subroutine


  !> Initializes the module. All arguments correspond to similarly
  !! named ones in the module description.
  subroutine particlewall_init2(alpha_a_, alpha_b_, LJdist_, &
       isuniformalignment_, sig_, alpha_lj_, sigwall_lj_, wall_density_, &
       epswall_lj_, eps_)
    real(dp), intent(in), optional :: alpha_a_, alpha_b_, LJdist_, sig_, &
         alpha_lj_, sigwall_lj_, wall_density_, epswall_lj_, eps_
    logical, intent(in), optional :: isuniformalignment_
    if (present(alpha_a_)) alpha_a = alpha_a_
    if (present(alpha_b_)) alpha_b = alpha_b_
    if (present(LJdist_)) LJdist = LJdist_
    if (present(isuniformalignment_)) isuniformalignment = isuniformalignment_
    if (present(sig_)) sig = sig_
    if (present(alpha_lj_)) alpha_lj = alpha_lj_
    if (present(sigwall_lj_)) sigwall_lj = sigwall_lj_
    if (present(wall_density_)) wall_density = wall_density_
    if (present(epswall_lj_)) epswall_lj = epswall_lj_
    if (present(eps_)) eps = eps_
  end subroutine


  !> Writes the parameters of this module with the format and output
  !! unit defined by @p writer.
  subroutine particlewall_writeparameters(writer)
    type(parameter_writer), intent(in) :: writer
    call writecomment(writer, 'General wall parameters')
    call writeparameter(writer, 'wall_density', wall_density)

    call writecomment(writer, 'GB-wall interaction parameters')
    call writeparameter(writer, 'alpha_A', alpha_a) 
    call writeparameter(writer, 'alpha_B', alpha_b) 
    call writeparameter(writer, 'LJ_dist', LJdist) 
    call writeparameter(writer, 'is_uniform_alignment', isuniformalignment)
    call writeparameter(writer, 'sigwall', sig)
    call writeparameter(writer, 'epswall', eps)

    call writecomment(writer, 'LJ-wall interaction parameters')
    call writeparameter(writer, 'alpha_LJ', alpha_lj)
    call writeparameter(writer, 'sigwall_LJ', sigwall_lj)
    call writeparameter(writer, 'epswall_LJ', epswall_lj)
  end subroutine


  !> Writes the parameters of this module with the format and output
  !! unit defined by @p writer.
  subroutine particlewall_to_json(json_val)
    type(json_value), intent(inout), pointer :: json_val
    
    call json_add(json_val, 'wall_density', wall_density)

    call json_add(json_val, 'rodwall', 'ljdimer-wall')
    call json_add(json_val, 'alpha_A', alpha_a) 
    call json_add(json_val, 'alpha_B', alpha_b) 
    call json_add(json_val, 'LJ_dist', LJdist) 
    call json_add(json_val, 'is_uniform_alignment', isuniformalignment)
    call json_add(json_val, 'sigwall', sig)
    call json_add(json_val, 'epswall', eps)

    call json_add(json_val, 'pointwall', 'ljcylinder')
    call json_add(json_val, 'alpha_LJ', alpha_lj)
    call json_add(json_val, 'sigwall_LJ', sigwall_lj)
    call json_add(json_val, 'epswall_LJ', epswall_lj)
  end subroutine


  !> Calculates the potential @p energy for a @p particle inside a
  !! cylindrical Lennard-Jones cavity. The radius of the cavity is
  !! defined by @p simbox. If the particle is too close to the wall or
  !! inside the wall @p overlap == true.
  pure subroutine particlewall_potential(particle, simbox, energy, err)
    type(particledat), intent(in) :: particle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    integer, intent(out) :: err
    logical :: overlap
    err = 0
    if (particle%rod) then
      call gbwall(particle, simbox, energy, overlap)
    else
      call ljwall(particle, simbox, energy, overlap)
   end if
   if (overlap) err = 1
  end subroutine


  !> Returns the force acting on @p particle by the wall of a
  !! cylindrical cavity. Radius of the cavity is defined by @p simbox.
  !! See module description for details.
  function particlewall_force(particle, simbox) result(f)
    type(particledat), intent(in) :: particle
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)
    if (particle%rod) then
      f = gbwall_force(particle, simbox)
    else
      f = ljwall_force(particle, simbox)
    end if
  end function


  !> Calculates the interaction energy of a Lennard-Jones (LJ) particle
  !! and the wall of a cylindrical cavity.
  !! 
  !! @param ljparticle the LJ particle.
  !! @param simbox the simulation box defining the radius of the
  !!        cavity.
  !! @param energy the interaction energy.
  !! @param overlap is true if the @p ljparticle has penetrated the
  !! wall too much.
  !! 
  pure subroutine ljwall(ljparticle, simbox, energy, overlap)
    implicit none
    type(particledat), intent(in) :: ljparticle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    real(dp) :: r
    real(dp) :: r_cylinder
    overlap = .false.
    energy = 0._dp
    r_cylinder = getx(simbox)/2._dp 
    r = sqrt(ljparticle%x**2 + ljparticle%y**2)
    if(r >= r_cylinder) then
      overlap = .true.
      return
    end if
    energy = ljcylinderpotential(epswall_lj, wall_density, sigwall_lj, alpha_lj, r, r_cylinder)
  end subroutine


  !> Returns the force exerted on a LJ particle @p ljparticle by the
  !! wall of the cylindrical cavity. The radius of the cavity is
  !! defined by the @p simbox.
  function ljwall_force(ljparticle, simbox) result(f)
    implicit none
    type(particledat), intent(in) :: ljparticle
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)

    real(dp) :: r
    real(dp) :: r_cylinder
    r_cylinder = getx(simbox)/2._dp 
    r = sqrt(ljparticle%x**2 + ljparticle%y**2)
    !! Set the direction of f:
    f = [ljparticle%x, ljparticle%y, 0._dp] / r
    !! What should be done when r is near zero in the division above?
    f = f * ljcylinderforce(epswall_lj, wall_density, sigwall_lj, alpha_lj, r, r_cylinder)
  end function


  !> Calculates the potential energy of a rodlike particle which
  !! interacts with the wall of a cylindrical cavity. The particle
  !! interacts with the wall via two embedded Lennard-Jones interaction
  !! sites.
  !! 
  !! @param gbparticle the rodlike particle.
  !! @param simbox the simulation box defining the dimensions of the
  !!        cavity.
  !! @param energy the interaction energy of the wall and the particle.
  !! @param overlap is true if the @p gbparticle has penetrated the
  !!        wall too much.
  !! 
  !! @see D. Micheletti et al. J. Chem. Phys. 123, 224705, 2005 for the
  !! interaction site model.
  !!
  pure subroutine gbwall(gbparticle, simbox, energy, overlap)
    implicit none
    type(particledat), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp), intent(out) :: energy
    logical, intent(out) :: overlap
    real(dp) :: rsite_a, rsite_b
    real(dp) :: fu, r_cylinder
    overlap = .false.
    energy = 0._dp
    r_cylinder = getx(simbox)/2._dp 
    call rarb(gbparticle, rsite_a, rsite_b)
    if(rsite_a >= r_cylinder .or. rsite_b >= r_cylinder) then
      overlap = .true.
      return
    end if
    if (isuniformalignment) then
      fu = angular(gbparticle)  
    else 
      fu = 1._dp
    end if
    energy = fu * (ljcylinderpotential(eps, wall_density, sig, alpha_a, rsite_a, r_cylinder) + &
    ljcylinderpotential(eps, wall_density, sig, alpha_b, rsite_b, r_cylinder))
  end subroutine gbwall


  !> Returns the force exerted on a GB particle @p gbparticle by the
  !! wall of a cylindrical cavity. The @p gbparticle interacts with the
  !! wall via two embedded LJ sites.
  function gbwall_force(gbparticle, simbox) result(f)
    type(particledat), intent(in) :: gbparticle
    type(poly_box), intent(in) :: simbox
    real(dp) :: f(3)
    real(dp) :: f_a(3), f_b(3)
    real(dp) :: rsite_a, rsite_b
    real(dp) :: fu, r_cylinder
    r_cylinder = getx(simbox) / 2. 
    call rarb(gbparticle, rsite_a, rsite_b)
    if (isuniformalignment) then
      fu = angular(gbparticle)  
    else 
      fu = 1.
    end if

    f_a = position(gbparticle) + orientation(gbparticle) * LJdist
    f_a(3) = 0._dp
    f_a = f_a / sqrt(dot_product(f_a, f_a))

    f_b = position(gbparticle) - orientation(gbparticle) * LJdist
    f_b(3) = 0._dp
    f_b = f_b / sqrt(dot_product(f_b, f_b))

    f = fu * (f_a * ljcylinderforce(eps, wall_density, sig, alpha_a, rsite_a, r_cylinder) + &
    f_b * ljcylinderforce(eps, wall_density, sig, alpha_b, rsite_b, r_cylinder))
  end function


  !> Returns the distances from cavity axis @p ra, @p rb for
  !! interaction sites embedded in @p particle.
  pure subroutine rarb(particle, ra, rb)
    type(particledat), intent(in) :: particle
    real(dp), intent(out) :: ra, rb
    real(dp) :: xa, ya, xb, yb
    xa = particle%x + LJdist * particle%ux
    ya = particle%y + LJdist * particle%uy
    xb = particle%x - LJdist * particle%ux
    yb = particle%y - LJdist * particle%uy
    ra = sqrt(xa**2 + ya**2)
    rb = sqrt(xb**2 + yb**2)
  end subroutine


  !> Returns the angular dependence of the potential with respect to the 
  !! cylinder axis (z-direction). 
  !!
  !! @param particle the particle to which the potential is calculated. 
  !! 
  pure real(dp) function angular(particle)
    type(particledat), intent(in) :: particle
    angular = (particle%uz)**2
  end function angular

end module particlewall
