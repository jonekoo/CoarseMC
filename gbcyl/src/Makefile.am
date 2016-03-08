##SUBDIRS=tests

## Programs to install
bin_PROGRAMS = \
  ptgbcyl \
  create_crystal \
  printpotentials


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
  set_inputparameters.py \
  update_json.py \
  convert_to_json.py \
  reorder_output.py \
  json_average.py

AM_FCFLAGS = \
-I${top_builddir}/mt_stream_f90-1.11 \
-I${top_builddir}/mt_stream_f90-1.11/f_jump_ahead_coeff \
-I${top_builddir}/claf90


utils_src = \
  num_kind.f90 \
  utils.f90 \
  m_fileunit.f90

model_src = \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  json_module.F90 \
  m_json_wrapper.f90 \
  process_not_found.inc \
  check_parameter.inc \
  class_poly_box.f90 \
  particle_mover.f90 \
  m_particle.f90 \
  m_point.f90 \
  m_rod.f90 



noinst_LTLIBRARIES = libmodel.la

libmodel_la_SOURCES = \
  $(utils_src) \
  $(model_src) \
  $(energy_src)

m_json_wrapper.$(OBJEXT) : m_json_wrapper.f90 check_parameter.inc process_not_found.inc

energy_src = \
  m_lj.f90 \
  m_gayberne.f90 \
  m_gblj.f90 \
  ljcylinder.f90 \
  particlewall.f90 \
  class_simplelist.f90

## Use dsyev from LAPACK
LDADD = $(ldadd_mtstream) $(ldadd_hyp2f1) $(top_builddir)/claf90/libclaf.la libmodel.la

## Use hyp2f1 from cephes
energy_src += m_hyp2f1_negint.f90 cylinder_integrals.f90 ljwall_cephes.f90
ldadd_hyp2f1 = $(top_builddir)/cephes/libhyp2f1.la

ldadd_mtstream = ${top_builddir}/mt_stream_f90-1.11/libmt_stream.la

AM_FCFLAGS += -I$(top_builddir)/cephes/cmath


ptgbcyl_SOURCES = \
  genvoltrial.f90 \
  beta_exchange.f90 \
  m_particlegroup.f90 \
  m_lj_interaction.f90 \
  m_gb_interaction.f90 \
  m_gblj_interaction.f90 \
  m_interaction_factory.f90 \
  m_particlejson_parser.f90 \
  m_nvt_engine.F90 \
  m_npt_engine.F90 \
  mc_engine.F90 \
  gbcyl.f90 

create_crystal_SOURCES = \
  crystal.f90 \
  parse_cmdline.inc \
  create_crystal.f90

printpotentials_SOURCES = \
  genvoltrial.f90 \
  beta_exchange.f90 \
  m_particlegroup.f90 \
  m_lj_interaction.f90 \
  m_gb_interaction.f90 \
  m_gblj_interaction.f90 \
  m_interaction_factory.f90 \
  m_particlejson_parser.f90 \
  m_nvt_engine.F90 \
  m_npt_engine.F90 \
  mc_engine.F90 \
  printpotentials.f90 

#printforces_SOURCES = \
#  printforces.f90 

dist_check_SCRIPTS = test_zerosweeps.py
TESTS = test_zerosweeps.py
TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_PY_LOG_FLAGS = -v
dist_check_DATA = tests/zero_sweeps/0/input-0.json tests/zero_sweeps/total_energies-reference.txt

test_zerosweeps.py : tests/zero_sweeps/0/restart-0.json tests/zero_sweeps/total_energies-reference.txt

tests/zero_sweeps/0/restart-0.json : tests/zero_sweeps/0/input-0.json $(top_builddir)/src/ptgbcyl
	mkdir -p tests/zero_sweeps/0;
	mpirun -np 1 $(top_builddir)/src/ptgbcyl -i $(top_srcdir)/src/tests/zero_sweeps/0/input-_I_.json -r tests/zero_sweeps/0/restart-_I_.json

if PFUNIT

driver_SOURCES = \
  class_parameterizer.f90 \
  class_parameter_writer.f90 \
  json_module.F90 \
  m_json_wrapper.f90 \
  m_gblj.f90 \
  m_gayberne.f90 \
  gayberne_pfunit.f90 \
  particle_mover.f90 \
  m_particle.f90 \
  class_poly_box.f90 \
  m_lj.f90 \
  class_simplelist.f90 \
  class_simplelist_pfunit.f90 \
  utils_pfunit.f90 \
  m_point.f90 \
  m_rod.f90 \
  m_particlejson_parser.f90 \
  particle_pfunit.f90 \
  driver.f90  


driver_FCFLAGS = -I${PFUNIT_DIR}/mod \
-I${top_builddir}/mt_stream_f90-1.11 \
-I${top_builddir}/mt_stream_f90-1.11/f_jump_ahead_coeff

driver_LDADD = ${PFUNIT_DIR}/lib/libpfunit.a ${ldadd_mtstream} libmodel.la

TESTS += driver 
check_PROGRAMS = driver

endif


MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)
	-rm -f $(bin_PROGRAMS)
	-rm -f tests/zero_sweeps/*/restart-0.json
	-rm -f tests/short/id-*/restart-0.json
	-rm -f tests/short/id-*/output-0.json
	-rm -f beta_exchange.stats

