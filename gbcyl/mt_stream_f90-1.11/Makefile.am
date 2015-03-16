JUMPLIBDIR = f_jump_ahead_coeff


noinst_LTLIBRARIES = libmt_stream.la
libmt_stream_la_SOURCES =mt_kind_defs.F90 ${JUMPLIBDIR}/gf2xe.F90 ${JUMPLIBDIR}/f_get_coeff.F90 mt_stream.F90


#check_stream_1 : check_stream_1.o mt_stream_debug.o $(JUMPOBJ)
#	$(FC) $(LDFLAGS) -static $^ -o $@ $(LIBS)

#check_stream_2 : check_stream_2.o mt_stream.o $(JUMPOBJ)
#	$(FC) $(LDFLAGS) -static $^ -o $@ $(LIBS)

#check_stream_3 : check_stream_3.o mt_stream_debug.o $(JUMPOBJ)
#	$(FC) $(LDFLAGS) -static $^ -o $@ $(LIBS)

#check_stream_1.o : check_stream.F90 mt_stream_debug.o
#	$(FC) $(FFLAGS) -D_DEBUG_ -D_CHECK_STREAM_ $(INCLUDES) -c $< -o $@

#check_stream_2.o : check_stream.F90 mt_stream.o
#	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

#check_stream_3.o : check_stream.F90 mt_stream_debug.o
#	$(FC) $(FFLAGS) -D_DEBUG_ -D_CHECK_SAVEREAD_ $(INCLUDES) -c $< -o $@

#mt_stream_debug.o : mt_stream.F90 mt_kind_defs.o
#	$(FC) $(FFLAGS) -D_DEBUG_ $(INCLUDES) -c $< -o $@

if MT_TESTS
check_PROGRAMS = check_stream_2
check_stream_2_LDADD = libmt_stream.la
check_stream_2_SOURCES = check_stream.F90

check_DATA = test_2.out sample_2.out
check_SCRIPTS = check_2.sh

test_2.out : check_stream_2
	./check_stream_2 > test_2.out

TESTS = check_2.sh
endif

MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)
