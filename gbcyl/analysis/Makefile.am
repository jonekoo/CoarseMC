## Programs to install
bin_PROGRAMS = \
  orientation_parameter \
  xe_nmr \
  xe_nmr_profile \
  xe_nmr_profile2 \
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
  cylindrical_layerorientation \
  print_gbxe_qcoupling \
  print_gbxe_shielding \
  print_xewall_coupling \
  print_xewall_shielding \
  print_xeslab_shielding \
  pair_distribution
#  virial_tensor 

MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT) $(bin_PROGRAMS)


ldadd_mtstream = ${top_builddir}/mt_stream_f90-1.11/libmt_stream.la

ldadd_hyp2f1 = $(top_builddir)/cephes/libhyp2f1.la

LDADD = $(top_builddir)/lapack/lapack_routine/libdsyev.la $(ldadd_mtstream) $(ldadd_hyp2f1) $(top_builddir)/claf90/libclaf.la $(top_builddir)/src/libmodel.la

AM_FCFLAGS = \
-I${top_builddir}/mt_stream_f90-1.11 \
-I${top_builddir}/mt_stream_f90-1.11/f_jump_ahead_coeff \
-I${top_builddir}/claf90
AM_FCFLAGS += -I$(top_builddir)/cephes/cmath
AM_FCFLAGS += -I$(top_builddir)/src

eigens_src = eigens_dsyev.f90
common_src =   m_particledat.f90 m_particle_factory.f90 

orientation_parameter_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  orientation_parameter.f90 

xe_nmr_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  xewall_shielding.f90 \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  lj_nmr.f90 \
  xe_nmr.f90 \
  ljcylinder_tensor.f90 

xe_nmr_profile_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  xewall_shielding.f90 \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  lj_nmr.f90 \
  xfunc_module.f90 \
  histogram.f90 \
  $(top_builddir)/src/parse_cmdline.inc \
  xe_nmr_profile.f90 \
  ljcylinder_tensor.f90 

xe_nmr_profile2_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  xewall_shielding.f90 \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  lj_nmr.f90 \
  xfunc_module.f90 \
  histogram.f90 \
  $(top_builddir)/src/parse_cmdline.inc \
  xe_nmr_profile2.f90 \
  ljcylinder_tensor.f90 

print_xewall_coupling_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  xewall_shielding.f90 \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  print_xewall_coupling.f90 \
  ljcylinder_tensor.f90 

print_xewall_shielding_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  xewall_shielding.f90 \
  m_constants.f90 \
  m_shielding.f90 \
  print_xewall_shielding.f90 \
  ljcylinder_tensor.f90 

print_xeslab_shielding_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  xewall_shielding.f90 \
  m_constants.f90 \
  m_shielding.f90 \
  print_xewall_shielding.f90 \
  ljslab_tensor.f90 

print_gbxe_qcoupling_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  xewall_shielding.f90 \
  m_constants.f90 \
  m_shielding.f90 \
  xewall_qcoupling.f90 \
  m_quadrupole_coupling.f90 \
  print_gbxe_qcoupling.f90 \
  ljcylinder_tensor.f90 

print_gbxe_shielding_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  m_rank2_tensor.f90 \
  xewall_shielding.f90 \
  m_constants.f90 \
  m_shielding.f90 \
  print_gbxe_shielding.f90 \
  ljcylinder_tensor.f90 

cylindrical_orientation_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  cylindrical_orientation.f90 \
  cylindrical_orientation_tensor.f90

cylindrical_layerorientation_SOURCES = \
  $(common_src) \
  orientational_ordering.f90 \
  ${eigens_src} \
  layernormal.f90 \
  cylindrical_layerorientation.f90 \
  cylindrical_layerorientationtensor.f90

povanimator_SOURCES = \
  $(common_src) \
  pov.f90 \
  povanimation.f90

tau1_SOURCES = \
  $(common_src) \
  tau1_negative.f90 \
  tau1_module.f90 \
  orientational_ordering.f90 \
  ${eigens_src} \
  tau1_program.f90 \
  simplemin.f90 \
  max_tau1.f90

pair_distribution_SOURCES = \
  $(common_src) \
  distribution.f90 \
  orientational_ordering.f90 \
  ${eigens_src} \
  pair_distribution.f90 

weighted_rdf_SOURCES = \
  $(common_src) \
  layernormal.f90 \
  ${eigens_src} \
  gr3dweighted.f90 \
  weighted_rdf.f90 

density_SOURCES = \
  $(common_src) \
  density.f90

psi6_SOURCES = \
  $(common_src) \
  psi6_module.f90 \
  orientational_ordering.f90 \
  ${eigens_src} \
  psi6.f90 

cylhist_SOURCES = \
  $(common_src) \
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
  $(common_src) \
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
  $(common_src) \
  simplemin.f90 \
  tau1_negative.f90 \
  tau1_module.f90 \
  layernormal.f90 \
  ${eigens_src} \
  smctau1.f90 \
  max_tau1.f90

smcpsi6_SOURCES = \
  $(common_src) \
  simplemin.f90 \
  psi6_module.f90 \
  layernormal.f90 \
  ${eigens_src} \
  smcpsi6.f90

#pressurezz_SOURCES = \
#  $(common_src) \
#  genvoltrial.f90 \
#  m_particlegroup.f90 \
#  beta_exchange.f90 \
#  mc_engine.f90 \
#  pressurezz.f90

#virial_tensor_SOURCES = \
#  $(common_src) \
#  class_pair_potential.f90 \
#  virial_tensor.f90

#onwallpressure_SOURCES = \
#  $(common_src) \
#  onwallpressure.f90


