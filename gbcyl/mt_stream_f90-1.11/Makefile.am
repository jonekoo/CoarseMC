JUMPLIBDIR = f_jump_ahead_coeff

#FCFLAGS := -Wall -fno-range-check -pedantic -Wconversion -Wunused-parameter -Warray-bounds -I/opt/local/include/openmpi -I$(JUMPLIBDIR) -ffast-math
#FCFLAGS := -I${JUMPLIBDIR}
FCFLAGS=-O3 -g -fno-range-check

noinst_LTLIBRARIES = libmt_stream.la
libmt_stream_la_SOURCES =mt_kind_defs.F90 ${JUMPLIBDIR}/gf2xe.F90 ${JUMPLIBDIR}/f_get_coeff.F90 mt_stream.F90

MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)
