noinst_LTLIBRARIES = libclaf.la

libclaf_la_SOURCES = kinds.f90 cla.f90

if CLAF_TESTS

check_PROGRAMS = cla_test

cla_test_SOURCES = cla_test.f90
cla_test_LDADD = libclaf.la
cla_test_AM_FCLAGS = -cpp

endif

MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)

