## Programs to install
bin_PROGRAMS = \
  ptgbcyl \
  printpotentials \
  printforces \
  orientation_parameter \
  cylindrical_orientation \
  tau1 \
  weighted_rdf \
  density \
  psi6 \
  povanimator \
  smctau1 \
  smcpsi6 \
  cylhist \
  cylhist2 \
  configurationconverter \
  cylindrical_layerorientation \
  pressurezz \
  virial_tensor \
  doublezsize \
  onwallpressure \
  remove_every_2nd \
  test_simplelist_with_openmp \
  pair_distribution \
  xe_nmr \
  xe_nmr_profile \
  xe_nmr_profile2 \
  test_gblj \
  print_gbxe_qcoupling \
  print_gbxe_shielding \
  print_xewall_coupling \
  print_xewall_shielding \
  print_xeslab_shielding \
  create_crystal \
  default_parameters



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
  slowcool.py \
  gbs2ljs.py \
  set_inputparameters.py



ldadd_mtstream = ${top_builddir}/mt_stream_f90-1.11/libmt_stream.la

AM_FCFLAGS = \
-I${top_builddir}/mt_stream_f90-1.11 \
-I${top_builddir}/mt_stream_f90-1.11/f_jump_ahead_coeff


model_src = \
  particle.f90 \
  class_poly_box.f90 \
  class_factory.f90

energy_src = \
  lj.f90 \
  gayberne.f90 \
  gblj.f90 \
  ljcylinder.f90 \
  particlewall.f90 \
  class_pair_potential.f90 \
  verlet.f90 \
  class_simplelist.f90 \
  energy.f90 

ptgbcyl_LDADD = ${ldadd_mtstream}

nr_src = num_kind.f90
## Use dsyev from LAPACK
LDADD = ${top_builddir}/lapack/lapack_routine/libdsyev.la
eigens_src = eigens_dsyev.f90
## Use hyp2f1 from cephes
energy_src += m_hyp2f1_negint.f90 cylinder_integrals.f90 ljwall_cephes.f90
ldadd_hyp2f1 = ${top_builddir}/cephes/libhyp2f1.la
ptgbcyl_LDADD += ${ldadd_hyp2f1}
printpotentials_LDADD = ${ldadd_hyp2f1}
printforces_LDADD = ${ldadd_hyp2f1}
AM_FCFLAGS += -I${top_builddir}/cephes/cmath
bin_PROGRAMS += test_hyp2f1 test_zhangI
test_hyp2f1_SOURCES =  m_hyp2f1_negint.f90 test_hyp2f1.f90
test_hyp2f1_LDADD = ${ldadd_hyp2f1}  ${ldadd_mtstream}
test_zhangI_SOURCES = utils.f90 m_hyp2f1_negint.f90 test_zhangI.f90 ljwall_cephes.f90
test_zhangI_LDADD = ${ldadd_hyp2f1}  ${ldadd_mtstream}
pressurezz_LDADD = ${ldadd_hyp2f1}
virial_tensor_LDADD = ${ldadd_hyp2f1}
onwallpressure_LDADD = ${ldadd_hyp2f1}
test_simplelist_with_openmp_LDADD = ${ldadd_hyp2f1}
xe_nmr_LDADD = ${LDADD} ${ldadd_hyp2f1}
xe_nmr_profile_LDADD = ${LDADD} ${ldadd_hyp2f1}
xe_nmr_profile2_LDADD = ${LDADD} ${ldadd_hyp2f1}
test_gblj_LDADD = ${ldadd_hyp2f1}
print_gbxe_qcoupling_LDADD = ${LDADD} ${ldadd_hyp2f1}
print_gbxe_shielding_LDADD = ${LDADD} ${ldadd_hyp2f1}
print_xewall_coupling_LDADD = ${LDADD} ${ldadd_hyp2f1}
print_xewall_shielding_LDADD = ${LDADD} ${ldadd_hyp2f1}
print_xeslab_shielding_LDADD = ${LDADD} ${ldadd_hyp2f1}

common_src = \
  ${nr_src} \
  rng.inc \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  particle_mover.f90 \
  genvoltrial.f90 

parallel_src = \
  mpi.f90 \
  beta_exchange.f90 \
  mc_sweep.f90 \
  mc_engine.f90 

## Below: Listings of source files needed by programs

ptgbcyl_SOURCES = \
  ${common_src} \
  ${energy_src} \
  ${parallel_src} \
  gbcyl.f90 

default_parameters_LDADD = $(ptgbcyl_LDADD)

default_parameters_SOURCES = \
  ${common_src} \
  ${energy_src} \
  ${parallel_src} \
  default_parameters.f90

create_crystal_SOURCES = \
  utils.f90 \
  ${model_src} \
  crystal.f90 \
  parse_cmdline.inc \
  create_crystal.f90

orientation_parameter_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  orientation_parameter.f90 

xe_nmr_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  cylinder_integrals.f90 \
  xewall_shielding.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  lj_nmr.f90 \
  xe_nmr.f90 \
  ljcylinder_tensor.f90 


xe_nmr_profile_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  cylinder_integrals.f90 \
  xewall_shielding.f90 \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  lj_nmr.f90 \
  xfunc_module.f90 \
  histogram.f90 \
  parse_cmdline.inc \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  xe_nmr_profile.f90 \
  ljcylinder_tensor.f90 

xe_nmr_profile2_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  cylinder_integrals.f90 \
  xewall_shielding.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  lj_nmr.f90 \
  xfunc_module.f90 \
  histogram.f90 \
  parse_cmdline.inc \
  xe_nmr_profile2.f90 \
  ljcylinder_tensor.f90 

print_xewall_coupling_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  cylinder_integrals.f90 \
  xewall_shielding.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  print_xewall_coupling.f90 \
  ljcylinder_tensor.f90 

print_xewall_shielding_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  cylinder_integrals.f90 \
  xewall_shielding.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  m_constants.f90 \
  m_shielding.f90 \
  print_xewall_shielding.f90 \
  ljcylinder_tensor.f90 

print_xeslab_shielding_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  cylinder_integrals.f90 \
  xewall_shielding.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  m_constants.f90 \
  m_shielding.f90 \
  print_xewall_shielding.f90 \
  ljslab_tensor.f90 

print_gbxe_qcoupling_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  cylinder_integrals.f90 \
  xewall_shielding.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  print_gbxe_qcoupling.f90 \
  ljcylinder_tensor.f90 

print_gbxe_shielding_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  cylinder_integrals.f90 \
  xewall_shielding.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  m_constants.f90 \
  m_shielding.f90 \
  print_gbxe_shielding.f90 \
  ljcylinder_tensor.f90 


remove_every_2nd_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  remove_every_2nd.f90

cylindrical_orientation_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  cylindrical_orientation.f90 \
  cylindrical_orientation_tensor.f90

cylindrical_layerorientation_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  orientational_ordering.f90 \
  ${eigens_src} \
  layernormal.f90 \
  cylindrical_layerorientation.f90 \
  cylindrical_layerorientationtensor.f90

povanimator_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  pov.f90 \
  povanimation.f90

tau1_SOURCES = \
  ${nr_src} \
  simplemin.f90 \
  m_fileunit.f90 \
  ${model_src} \
  tau1_negative.f90 \
  utils.f90 \
  tau1_module.f90 \
  orientational_ordering.f90 \
  ${eigens_src} \
  tau1_program.f90 \
  max_tau1.f90

pair_distribution_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  distribution.f90 \
  orientational_ordering.f90 \
  ${eigens_src} \
  pair_distribution.f90 

weighted_rdf_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  layernormal.f90 \
  ${eigens_src} \
  gr3dweighted.f90 \
  weighted_rdf.f90 

density_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  density.f90

psi6_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  psi6_module.f90 \
  orientational_ordering.f90 \
  ${eigens_src} \
  psi6.f90 

cylhist_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  xfunc_module.f90 \
  histogram.f90 \
  orientational_ordering.f90 \
  ${eigens_src} \
  simplemin.f90 \
  tau1_negative.f90 \
  tau1_module.f90 \
  psi6_module.f90 \
  layernormal.f90 \
  distribution.f90 \
  parse_cmdline.inc \
  cylhist.f90 \
  max_tau1.f90

cylhist2_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  xfunc_module.f90 \
  histogram.f90 \
  orientational_ordering.f90 \
  ${eigens_src} \
  simplemin.f90 \
  tau1_negative.f90 \
  tau1_module.f90 \
  psi6_module.f90 \
  layernormal.f90 \
  distribution.f90 \
  parse_cmdline.inc \
  cylhist2.f90 \
  max_tau1.f90  

smctau1_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  simplemin.f90 \
  tau1_negative.f90 \
  tau1_module.f90 \
  layernormal.f90 \
  ${eigens_src} \
  smctau1.f90 \
  max_tau1.f90

smcpsi6_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  simplemin.f90 \
  psi6_module.f90 \
  layernormal.f90 \
  ${eigens_src} \
  smcpsi6.f90

printpotentials_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  printpotentials.f90 

printforces_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  ${model_src} \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  printforces.f90 

configurationconverter_SOURCES = \
  ${nr_src} \
  ${model_src} \
  utils.f90 \
  m_fileunit.f90 \
  configurationconverter.f90

test_gblj_SOURCES = \
  ${nr_src} \
  ${model_src} \
  utils.f90 \
  m_fileunit.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  test_gblj.f90 

doublezsize_SOURCES = \
  ${nr_src} \
  ${model_src} \
  utils.f90 \
  m_fileunit.f90 \
  doublezsize.f90

pressurezz_SOURCES = \
  ${nr_src} \
  utils.f90 \
  m_fileunit.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${model_src} \
  ${energy_src} \
  pressurezz.f90

virial_tensor_SOURCES = \
  ${nr_src} \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${model_src} \
  ${energy_src} \
  utils.f90 \
  m_fileunit.f90 \
  virial_tensor.f90

onwallpressure_SOURCES = \
  ${nr_src} \
  m_fileunit.f90 \
  ${model_src} \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  utils.f90 \
  onwallpressure.f90

test_simplelist_with_openmp_SOURCES = \
  ${nr_src} \
  ${model_src} \
  m_fileunit.f90 \
  utils.f90 \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  ${energy_src} \
  test_simplelist_with_openmp.f90

testrunner_SOURCES = \
  num_kind.f90 \
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
  num_kind.f90 \
  m_fileunit.f90 \
  class_parameterizer.f90 \
  rng.inc mtfort90.f90 \
  utils.f90 \
  class_parameter_writer.f90 \
  gblj.f90 \
  gayberne.f90 \
  gayberne_pfunit.f90 \
  particle.f90 \
  class_poly_box.f90 \
  lj.f90 \
  class_pair_potential.f90 \
  class_simplelist.f90 \
  class_simplelist_pfunit.f90 \
  driver.F90  

PFUNIT_DIR=${top_builddir}/pFUnit

driver_FCFLAGS = ${FCFLAGS} -I${PFUNIT_DIR}/mod \
	-I${top_builddir}/mt_stream_f90-1.11
driver_LDADD = ${PFUNIT_DIR}/lib/libpfunit.a ${ldadd_mtstream}

testrunner_FCFLAGS = ${FCFLAGS} -I${top_builddir}/ftnunit \
	-I${top_builddir}/mt_stream_f90-1.11

testrunner_LDADD = ${top_builddir}/ftnunit/libftnunit.a ${ldadd_mtstream}

TESTS = driver
dist_check_SCRIPTS = ${top_builddir}/ftnunit/runtests.sh
check_PROGRAMS = testrunner driver 

MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)
	-rm -f *_ser.f90

