dist_check_SCRIPTS = test_short.py
TESTS = test_short.py
TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_PY_LOG_FLAGS = -v
dist_check_DATA = shortreference-0.json shortreference-1.json shortinput-0.json shortinput-1.json

test_short.py : shortrestart-0.json shortrestart-1.json shortreference-0.json shortreference-1.json 

shortrestart-0.json shortrestart-1.json : shortinput-0.json shortinput-1.json $(top_builddir)/src/ptgbcyl
	mpirun -np 2 $(top_builddir)/src/ptgbcyl -i $(top_srcdir)/tests/shortinput-_I_.json -o shortoutput-_I_.json -r shortrestart-_I_.json

mostlyclean-local :
	-rm -f shortrestart-0.json shortrestart-1.json shortoutput-0.json shortoutput-1.json