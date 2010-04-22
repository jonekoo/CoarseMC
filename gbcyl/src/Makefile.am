bin_PROGRAMS = ptgbcyl printpotentials #tracer #gbcyl

printpotentials_SOURCES = \
  nrtype.f90 \
  mtfort90.f90 \
  utils.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  gayberne.f90 \
  particle.f90 \
  cylinder.f90 \
  box.f90 \
  class_poly_box.f90 \
  particlewall.f90 \
  horner.f90 \
  nrutil.f90 \
  rf_s.f90 \
  rd_s.f90 \
  printpotentials.f90

#tracer_SOURCES = \
  ${common_src} \
  cell_energy.f90 \
  class_geometryreader.f90 \
  class_particlesystem.f90 \
  class_nptintegral.f90 \
  tracer.f90

common_src = \
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
  class_cylformatter.f90 \
  gayberne.f90 \
  cell.f90 \
  particlewall.f90 \
  class_pair_potential.f90 \
  verlet.f90 \
  all_pairs.f90 \
  class_nbrcelliterator.f90 \
  cell_energy.f90 \
  class_poly_nbrlist.f90 \
  energy.f90 \
  horner.f90 \
  rd_s.f90 \
  rd_v.f90 \
  rf_s.f90 \
  rf_v.f90 

mpi_src = \
  mpi.f90 \
  io.f90 \
  pt.f90 \
  mc_sweep.f90 \
  mc_engine.f90 

##%_ser.f90 : %.f90
##	awk '!(/mpi/ || /MPI/ || /pt_/) {print}' $< > $@

ptgbcyl_SOURCES = \
  ${common_src} \
  ${mpi_src} \
  gbcyl.f90

TESTS_ENVIRONMENT = openmpirun -np 2  
TESTS = TestRunner
nodist_TESTS = TestRunner parRestartTest1.sh parRestartTest2.sh
#dist_check_SCRIPTS = initTest.sh parRestartTest1.sh parRestartTest2.sh cleanupTest.sh
#dist_data_DATA = testconfiguration.txt
#dist_check_DATA = testconfiguration.txt   ../data/ptgbcyl.in ../data/r9-n1290-cr.gb
check_PROGRAMS = TestRunner
TestRunner_SOURCES = \
  ${common_src} \
  ${mpi_src} \
  lj.f90 \
  configurations.f90 \
  testutils.f90 \
  class_parameterizer_fun.f90 \
  verlet_fun.f90 \
  cell_energy_fun.f90 \
  class_cylformatter_fun.f90 \
  particle_fun.f90 \
  box_fun.f90 \
  mtmod_fun.f90 \
  gayberne_fun.f90 \
  io_fun.f90 \
  cell_fun.f90 \
  cylinder_fun.f90 \
  pt_fun.f90 \
  class_nbrcelliterator_fun.f90 \
  class_poly_nbrlist_fun.f90 \
  TestRunner.f90

MODEXT = @MODEXT@


mostlyclean-local :
	-rm -f *.$(MODEXT)
	-rm -f *_ser.f90
