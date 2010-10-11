bin_PROGRAMS = \
  ptgbcyl \
  printpotentials \
  createbox \
  cutcylinder \
  smctau1 \
  orientation_parameter \
  cylindrical_orientation \
  tau1 \
  radial_distribution \
  density \
  psi6 \
  povmaker \
  particlehistogram \
  povanimator \
  smctau1

dist_bin_SCRIPTS = \
  rdf_average.sh \
  loop-script-over-args.sh \
  remove-non-numeric-lines.sh \
  density_cyl.awk \
  show-non-numeric-lines.sh \
  select.awk \
  periodic_select.awk \
  volume_cyl.awk \
  find_closest.awk \
  dots2commas.sh \
  every-xth-field.awk \
  pipe-columns.awk \
  pt.sh \
  createinput.sh

orientation_parameter_SOURCES = \
  ${nr_src} \
  mtfort90.f90 \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  state_reader.f90 \
  orientational_ordering.f90 \
  orientation_parameter_wrapper.f90 \
  orientation_parameter.f90

cylindrical_orientation_SOURCES = \
  ${nr_src} \
  utils.f90 \
  state_reader.f90 \
  m_fileunit.f90 \
  mtfort90.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  cylindrical_orientation.f90 \
  cylindrical_orientation_wrapper.f90 \
  cylindrical_orientation_tensor.f90

povanimator_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  mtfort90.f90 \
  ${model_src} \
  state_reader.f90 \
  pov.f90 \
  povanimation.f90

tau1_SOURCES = \
  ${nr_src} \
  simplemin.f90 \
  m_fileunit.f90 \
  mtfort90.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  tau1_negative.f90 \
  utils.f90 \
  tau1_module.f90 \
  state_reader.f90 \
  tau1_program.f90 

radial_distribution_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  mtfort90.f90 \
  ${model_src} \
  state_reader.f90 \
  radial_distribution_wrapper.f90 \
  radial_distribution.f90 \
  gr3d.f90

density_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  mtfort90.f90 \
  ${model_src} \
  state_reader.f90 \
  density.f90

psi6_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  mtfort90.f90 \
  ${model_src} \
  state_reader.f90 \
  psi6_module.f90 \
  orientational_ordering.f90 \
  psi6.f90

povmaker_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  mtfort90.f90 \
  ${model_src} \
  io.f90 \
  state_reader.f90 \
  povize.f90 
                   
particlehistogram_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  mtfort90.f90 \
  ${model_src} \
  state_reader.f90 \
  histogram.f90 \
  orientational_ordering.f90 \
  simplemin.f90 \
  tau1_negative.f90 \
  tau1_module.f90 \
  xfunc_module.f90 \
  particlehistogram.f90 
                            
nr_src = \
  nrtype.f90 \
  nrutil.f90 \
  brent.f90 \
  rf_s.f90 \
  rf_v.f90 \
  rd_s.f90 \
  rd_v.f90 \
  nr.f90 

smctau1_SOURCES = \
  ${nr_src} \
  mtfort90.f90 \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  simplemin.f90 \
  orientational_ordering.f90 \
  tau1_negative.f90 \
  tau1_module.f90 \
  layernormal.f90 \
  smctau1.f90

printpotentials_SOURCES = \
  ${nr_src} \
  mtfort90.f90 \
  utils.f90 \
  m_fileunit.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  gayberne.f90 \
  particle.f90 \
  cylinder.f90 \
  box.f90 \
  class_poly_box.f90 \
  particlewall.f90 \
  horner.f90 \
  printpotentials.f90

createbox_SOURCES = \
  nrtype.f90 \
  utils.f90 \
  m_fileunit.f90 \
  class_parameter_writer.f90 \
  class_parameterizer.f90 \
  mtfort90.f90 \
  particle.f90 \
  createbox.f90 

cutcylinder_SOURCES = \
  ${nr_src} \
  io.f90 \
  utils.f90 \
  mtfort90.f90 \
  m_fileunit.f90 \
  ${model_src} \
  tightpack.f90 

#tracer_SOURCES = \
  ${common_src} \
  cell_energy.f90 \
  class_geometryreader.f90 \
  class_particlesystem.f90 \
  class_nptintegral.f90 \
  tracer.f90

model_src = \
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
  horner.f90 

common_src = \
  ${nr_src} \
  mtfort90.f90 \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} 

mpi_src = \
  mpi.f90 \
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
  cell_fun.f90 \
  cylinder_fun.f90 \
  pt_fun.f90 \
  class_nbrcelliterator_fun.f90 \
  class_poly_nbrlist_fun.f90 \
  class_poly_box_fun.f90 \
  TestRunner.f90

MODEXT = @MODEXT@


mostlyclean-local :
	-rm -f *.$(MODEXT)
	-rm -f *_ser.f90

