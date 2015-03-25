program test_particlegroup
use class_gblj
use m_rodgroup
use m_spheregroup
use m_gblj

!! create particlegroups
class(rodgroup), allocatable, target :: solvent
class(spheregroup), allocatable, target :: solute

!! connect the groups with interactions
call solvent%add_interaction(gblj_rods(solute))
call solute%add_interaction(gblj_spheres(solvent))
!! run simulation

end program
