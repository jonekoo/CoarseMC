dist_check_SCRIPTS = test_short.py test_zerosweeps.py
TESTS = test_short.py test_zerosweeps.py
TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_PY_LOG_FLAGS = -v
dist_check_DATA = shortreference-0.json shortreference-1.json shortinput-0.json shortinput-1.json total_energies-reference.txt zerosweepsinput-0.json

test_short.py : shortrestart-0.json shortrestart-1.json shortreference-0.json shortreference-1.json 

shortrestart-0.json shortrestart-1.json : shortinput-0.json shortinput-1.json $(top_builddir)/src/ptgbcyl
	mpirun -np 2 $(top_builddir)/src/ptgbcyl -i $(top_srcdir)/tests/shortinput-_I_.json -o shortoutput-_I_.json -r shortrestart-_I_.json

test_zerosweeps.py : zerosweepsrestart-0.json total_energies-reference.txt

zerosweepsrestart-0.json : zerosweepsinput-0.json $(top_builddir)/src/ptgbcyl
	mpirun -np 1 $(top_builddir)/src/ptgbcyl -i $(top_srcdir)/tests/zerosweepsinput-_I_.json -o zerosweepsoutput-_I_.json -r zerosweepsrestart-_I_.json


mostlyclean-local :
	-rm -f shortrestart-0.json shortrestart-1.json shortoutput-0.json shortoutput-1.json
	-rm -f zerosweepsrestart-0.json zerosweepsoutput-0.json

