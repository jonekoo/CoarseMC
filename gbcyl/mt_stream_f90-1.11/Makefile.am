JUMPLIBDIR = f_jump_ahead_coeff

noinst_LTLIBRARIES = libmt_stream.la
libmt_stream_la_SOURCES =mt_kind_defs.F90 ${JUMPLIBDIR}/gf2xe.F90 ${JUMPLIBDIR}/f_get_coeff.F90 mt_stream.F90

MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)
