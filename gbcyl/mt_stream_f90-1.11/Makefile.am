JUMPLIBDIR = f_jump_ahead_coeff


noinst_LTLIBRARIES = libmt_stream.la
libmt_stream_la_SOURCES =mt_kind_defs.F90 ${JUMPLIBDIR}/gf2xe.F90 ${JUMPLIBDIR}/f_get_coeff.F90 mt_stream.F90

TEST_EXTENSIONS = .sh

if MT_TESTS

TESTS = check_1.sh check_2.sh check_3.sh

SH_LOG_COMPILER = $(BASH)

check_PROGRAMS = check_stream_1 check_stream_2 check_stream_3
check_DATA = sample_1.out test_2.out sample_2.out
check_SCRIPTS = check_1.sh check_2.sh check_3.sh

check_LTLIBRARIES = libmt_stream_debug.la
libmt_stream_debug_la_SOURCES = $(libmt_stream_la_SOURCES)
libmt_stream_debug_la_FCFLAGS = -D_DEBUG_

check_stream_1_SOURCES = check_stream.F90
check_stream_1_FCFLAGS = -D_DEBUG_ -D_CHECK_STREAM_
check_stream_1_LDADD = libmt_stream_debug.la

check_stream_2_SOURCES = check_stream.F90
check_stream_2_LDADD = libmt_stream.la

check_stream_3_SOURCES = check_stream.F90
check_stream_3_FCFLAGS = -D_DEBUG_ -D_CHECK_SAVEREAD_
check_stream_3_LDADD = libmt_stream_debug.la


test_2.out : check_stream_2
	./check_stream_2 > test_2.out

endif

MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)

