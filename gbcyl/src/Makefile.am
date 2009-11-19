bin_PROGRAMS = ptgbcyl 
ptgbcyl_SOURCES = \
  nrtype.f90 \
  nrutil.f90 \
  nr.f90 \
  mtfort90.f90 \
  utils.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  particle.f90 \
  box.f90 \
  cylinder.f90 \
  class_poly_box.f90 \
  gayberne.f90 \
  cell.f90 \
  particlewall.f90 \
  verlet.f90 \
  all_pairs.f90 \
  energy.f90 \
  horner.f90 \
  mpi.f90 \
  io.f90 \
  pt.f90 \
  mc_sweep.f90 \
  rd_s.f90 \
  rd_v.f90 \
  rf_s.f90 \
  rf_v.f90 \
  mc_engine.f90 \
  gbcyl.f90

TESTS = TestRunner
dist_check_SCRIPTS = restartTest.sh
check_PROGRAMS = TestRunner
TestRunner_SOURCES = nrtype.f90 \
                     io.f90 \
                     io_fun.f90 \
                     class_parameter_writer.f90 \
                     class_parameterizer.f90 \
                     class_parameterizer_fun.f90 \
                     gayberne.f90 \
                     gayberne_fun.f90 \
                     mpi.f90 \
                     mtfort90.f90 \
                     utils.f90 \
                     cylinder.f90 \
                     particle.f90 \
                     particle_fun.f90 \
                     box.f90 \
                     box_fun.f90 \
                     class_poly_box.f90 \
                     pt.f90 \
                     configurations.f90 \
                     pt_fun.f90 \
                     lj.f90 \
                     verlet.f90 \
                     verlet_fun.f90 \
                     cell.f90 \
                     cell_fun.f90 \
                     TestRunner.f90
MODEXT = @MODEXT@


mostlyclean-local :
	-rm -f *.$(MODEXT)
