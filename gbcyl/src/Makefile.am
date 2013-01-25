## Programs to install
bin_PROGRAMS = \
  ptgbcyl \
  printpotentials \
  createbox \
  smctau1 \
  orientation_parameter \
  cylindrical_orientation \
  tau1 \
  radial_distribution \
  weighted_rdf \
  density \
  psi6 \
  povanimator \
  smctau1 \
  smcpsi6 \
  cylhist \
  printgbgb \
  printljcylinder \
  configurationconverter \
  cylindrical_layerorientation \
  pressurezz \
  virialpressure \
  doublezsize \
  onwallpressure \
  1Ddf \
  remove_every_2nd

## Scripts to install
dist_bin_SCRIPTS = \
  rdf_average.sh \
  loop-script-over-args.sh \
  remove-non-numeric-lines.sh \
  show-non-numeric-lines.sh \
  find_closest.awk \
  pipe-columns.awk \
  makerestart.sh \
  replaceparameter.py \
  slowcool.py

##ftnunit : 
##	$(MAKE) -C ftnunit

##LIBS = 

## Listings of source files needed by programs

ldadd_mtstream = ${top_builddir}/mt_stream_f90-1.11/libmt_stream.la

AM_FCFLAGS = -I${top_builddir}/mt_stream_f90-1.11 -I${top_builddir}/mt_stream_f90-1.11/f_jump_ahead_coeff

nr_src = nrtype.f90 nrutil.f90 nr.f90 rf_s.f90 rd_s.f90 brent.f90 pythag.f90 tqli.f90 tred2.f90 jacobi.f90

orientation_parameter_SOURCES = \
  ${nr_src} \
  rng.inc mtfort90.f90  \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  class_factory.f90 \
  jacobiwrapper.f90 \
  orientational_ordering.f90 \
  orientation_parameter.f90

remove_every_2nd_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  class_factory.f90 \
  remove_every_2nd.f90

remove_every_2nd_LDADD = $(ldadd_mtstream)

cylindrical_orientation_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  jacobiwrapper.f90 \
  orientational_ordering.f90 \
  cylindrical_orientation.f90 \
  class_factory.f90 \
  cylindrical_orientation_tensor.f90

cylindrical_layerorientation_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  jacobiwrapper.f90 \
  orientational_ordering.f90 \
  layernormal.f90 \
  cylindrical_layerorientation.f90 \
  class_factory.f90 \
  cylindrical_layerorientationtensor.f90

povanimator_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  class_factory.f90 \
  pov.f90 \
  povanimation.f90

tau1_SOURCES = \
  ${nr_src} \
  simplemin.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  tau1_negative.f90 \
  utils.f90 \
  tau1_module.f90 \
  jacobiwrapper.f90 \
  orientational_ordering.f90 \
  class_factory.f90 \
  tau1_program.f90 

radial_distribution_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  class_factory.f90 \
  radial_distribution.f90 \
  gr3d.f90

1Ddf_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  class_factory.f90 \
  distribution.f90 \
  1Ddf.f90

weighted_rdf_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  class_factory.f90 \
  layernormal.f90 \
  gr3dweighted.f90 \
  weighted_rdf.f90 

density_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  class_factory.f90 \
  density.f90

psi6_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  class_factory.f90 \
  psi6_module.f90 \
  jacobiwrapper.f90 \
  orientational_ordering.f90 \
  psi6.f90 

cylhist_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  rng.inc mtfort90.f90  \
  ${model_src} \
  state_reader.f90 \
  histogram.f90 \
  jacobiwrapper.f90 \
  orientational_ordering.f90 \
  simplemin.f90 \
  tau1_negative.f90 \
  tau1_module.f90 \
  xfunc_module.f90 \
  psi6_module.f90 \
  layernormal.f90 \
  class_factory.f90 \
  bin_indices.f90 \
  distribution.f90 \
  cylhist.f90 

smctau1_SOURCES = \
  ${nr_src} \
  rng.inc mtfort90.f90  \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  simplemin.f90 \
  tau1_negative.f90 \
  tau1_module.f90 \
  layernormal.f90 \
  class_factory.f90 \
  smctau1.f90

smcpsi6_SOURCES = \
  ${nr_src} \
  rng.inc mtfort90.f90  \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  simplemin.f90 \
  psi6_module.f90 \
  layernormal.f90 \
  class_factory.f90 \
  smcpsi6.f90

printpotentials_SOURCES = \
  ${nr_src} \
  rng.inc mtfort90.f90  \
  utils.f90 \
  m_fileunit.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  gayberne.f90 \
  ljcylinder.f90 \
  horner.f90 \
  printpotentials.f90 

printpotentials_LDADD = ${ldadd_mtstream}
createbox_LDADD = ${ldadd_mtstream}
smctau1_LDADD = ${ldadd_mtstream}
orientation_parameter_LDADD = ${ldadd_mtstream}
cylindrical_orientation_LDADD = ${ldadd_mtstream}
tau1_LDADD = ${ldadd_mtstream}
radial_distribution_LDADD = ${ldadd_mtstream}
weighted_rdf_LDADD = ${ldadd_mtstream}
density_LDADD = ${ldadd_mtstream}
psi6_LDADD = ${ldadd_mtstream}
povanimator_LDADD = ${ldadd_mtstream}
smcpsi6_LDADD = ${ldadd_mtstream}
cylhist_LDADD = ${ldadd_mtstream}
printgbgb_LDADD = ${ldadd_mtstream}
printljcylinder_LDADD = ${ldadd_mtstream}
configurationconverter_LDADD = ${ldadd_mtstream}
cylindrical_layerorientation_LDADD = ${ldadd_mtstream}
1Ddf_LDADD=${ldadd_mtstream}

printljcylinder_SOURCES = \
  ${nr_src} \
  rng.inc mtfort90.f90  \
  utils.f90 \
  ljcylinder.f90 \
  horner.f90 \
  printljcylinder.f90 

printgbgb_SOURCES = \
  nrtype.f90 \
  rng.inc mtfort90.f90  \
  utils.f90 \
  m_fileunit.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  gayberne.f90 \
  printgbgb.f90 

createbox_SOURCES = \
  nrtype.f90 \
  utils.f90 \
  m_fileunit.f90 \
  class_parameter_writer.f90 \
  class_parameterizer.f90 \
  rng.inc mtfort90.f90  \
  particle.f90 \
  class_poly_box.f90 \
  class_factory.f90 \
  createbox.f90 

configurationconverter_SOURCES = \
  ${nr_src} \
  ${model_src} \
  utils.f90 \
  rng.inc mtfort90.f90  \
  m_fileunit.f90 \
  class_factory.f90 \
  configurationconverter.f90

doublezsize_SOURCES = \
  ${nr_src} \
  ${model_src} \
  utils.f90 \
  rng.inc mtfort90.f90  \
  m_fileunit.f90 \
  class_factory.f90 \
  doublezsize.f90

#tracer_SOURCES = \
  ${common_src} \
  class_geometryreader.f90 \
  class_particlesystem.f90 \
  class_nptintegral.f90 \
  tracer.f90

model_src = \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  particle.f90 \
  class_poly_box.f90 \
  lj.f90 \
  gayberne.f90 \
  ljcylinder.f90 \
  particlewall.f90 \
  class_pair_potential.f90 \
  verlet.f90 \
  all_pairs.f90 \
  cell.f90 \
  class_simplelist.f90 \
  class_poly_nbrlist.f90 \
  energy.f90 \
  genvoltrial.f90 \
  horner.f90 

common_src = \
  ${nr_src} \
  rng.inc \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} 

parallel_src = \
  mpi.f90 \
  pt.f90 \
  mc_sweep.f90 \
  mc_engine.f90 

pressurezz_SOURCES = \
  ${common_src} \
  class_factory.f90 \
  ${parallel_src} \
  pressurezz.f90

pressurezz_LDADD = ${ldadd_mtstream}

doublezsize_LDADD = ${ldadd_mtstream}

onwallpressure_LDADD=${ldadd_mtstream}

virialpressure_SOURCES = \
  ${common_src} \
  class_factory.f90 \
  ${parallel_src} \
  virialpressure.f90

onwallpressure_SOURCES = \
  ${common_src} \
  class_factory.f90 \
  onwallpressure.f90

virialpressure_LDADD = ${ldadd_mtstream}

ptgbcyl_SOURCES = \
  ${common_src} \
  class_factory.f90 \
  ${parallel_src} \
  gbcyl.f90 

ptgbcyl_LDADD = ${ldadd_mtstream}

testrunner_SOURCES = \
  nrtype.f90 \
  m_fileunit.f90 \
  utils.f90 \
  utils_test.f90 \
  class_parameterizer.f90 \
  rng.inc mtfort90.f90  \
  class_parameter_writer.f90 \
  gayberne.f90 \
  gayberne_test.f90 \
  particle.f90 \
  particle_test.f90 \
  class_poly_box.f90 \
  class_poly_box_test.f90 \
  class_factory.f90 \
  class_factory_test.f90 \
  testrunner.f90 

driver_SOURCES = \
  Vectoroperations.F90 \
  nrtype.f90 \
  m_fileunit.f90 \
  class_parameterizer.f90 \
  rng.inc mtfort90.f90 \
  utils.f90 \
  class_parameter_writer.f90 \
  gayberne.f90 \
  gayberne_pfunit.f90 \
  particle.f90 \
  class_poly_box.f90 \
  cell.f90\
  class_simplelist.f90 \
  class_simplelist_pfunit.f90 \
  tests.F90 \
  driver.F90  

PFUNIT_DIR=/Users/jounikarjalainen/pFUnit

driver_FCFLAGS = ${FCFLAGS} -I${PFUNIT_DIR}/mod  -I${top_builddir}/mt_stream_f90-1.11
driver_LDADD = ${PFUNIT_DIR}/lib/libpfunit.a ${ldadd_mtstream}

testrunner_FCFLAGS = ${FCFLAGS} -I${top_builddir}/ftnunit -I${top_builddir}/mt_stream_f90-1.11
testrunner_LDADD = ${top_builddir}/ftnunit/libftnunit.a ${ldadd_mtstream}

#ftnunit_TESTS=testrunner
#ftnunit_TESTS_ENVIRONMENT=sh ${top_srcdir}/ftnunit/runtests.sh

##TESTS_ENVIRONMENT = openmpirun -np 2  
#TESTS_ENVIRONMENT = sh ${top_srcdir}/ftnunit/runtests.sh
#TESTS_ENVIRONMENT = mpirun -np 2 
TESTS = driver # testrunner 
##TESTS = driver # testrunner #test_gayberne test_particle
##nodist_TESTS = TestRunner parRestartTest1.sh parRestartTest2.sh
#dist_check_SCRIPTS = initTest.sh parRestartTest1.sh parRestartTest2.sh cleanupTest.sh
dist_check_SCRIPTS = ${top_builddir}/ftnunit/runtests.sh
#dist_data_DATA = testconfiguration.txt
#dist_check_DATA = testconfiguration.txt   ../data/ptgbcyl.in ../data/r9-n1290-cr.gb
##check_PROGRAMS = TestRunner
check_PROGRAMS = testrunner driver #test_gayberne test_particle #TestRunner

##distclean-local :
##	rm runtests.log

#TestRunner_SOURCES = \
#  ${common_src} \
#  ${mpi_src} \
#  lj.f90 \
#  configurations.f90 \
#  testutils.f90 \
#  class_parameterizer_fun.f90 \
#  verlet_fun.f90 \
#  class_cylformatter_fun.f90 \
#  particle_fun.f90 \
#  box_fun.f90 \
#  mtmod_fun.f90 \
#  gayberne_fun.f90 \
#  cylinder_fun.f90 \
#  pt_fun.f90 \
#  class_nbrcelliterator_fun.f90 \
#  class_poly_nbrlist_fun.f90 \
#  class_poly_box_fun.f90 \
#  TestRunner.f90

MODEXT = @MODEXT@


mostlyclean-local :
	-rm -f *.$(MODEXT)
	-rm -f *_ser.f90

