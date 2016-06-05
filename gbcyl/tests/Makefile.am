dist_check_SCRIPTS = test_short.py test_zerosweeps.py test_samename.py test_missinginput.py
TESTS = test_short.py test_zerosweeps.py test_samename.py test_missinginput.py
TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_PY_LOG_FLAGS = -v
dist_check_DATA = shortreference-0.json shortreference-1.json shortinput-0.json shortinput-1.json total_energies-reference.txt zerosweepsinput-0.json samename-0.json

export RUN_CMD=mpirun -np

input-0.json : $(top_srcdir)/tests/samename-0.json
	cp $(top_srcdir)/tests/samename-0.json input-0.json

test_samename.py : input-0.json

test_short.py : shortrestart-0.json shortrestart-1.json shortreference-0.json shortreference-1.json 

test_missinginput.json : $(top_srcdir)/tests/shortinput-0.json
	cp $(top_srcdir)/tests/shortinput-0.json test_missinginput.json

test_missinginput.py : test_missinginput.json

shortrestart-0.json shortrestart-1.json : shortinput-0.json shortinput-1.json $(top_builddir)/src/ptgbcyl
	$(RUN_CMD) 2 $(top_builddir)/src/ptgbcyl -i $(top_srcdir)/tests/shortinput-_I_.json -o shortoutput-_I_.json -r shortrestart-_I_.json

test_zerosweeps.py : zerosweepsrestart-0.json total_energies-reference.txt

zerosweepsrestart-0.json : zerosweepsinput-0.json $(top_builddir)/src/ptgbcyl
	$(RUN_CMD) 1 $(top_builddir)/src/ptgbcyl -i $(top_srcdir)/tests/zerosweepsinput-_I_.json -o zerosweepsoutput-_I_.json -r zerosweepsrestart-_I_.json

mostlyclean-local :
	-rm -f shortrestart-0.json shortrestart-1.json shortoutput-0.json shortoutput-1.json
	-rm -f zerosweepsrestart-0.json zerosweepsoutput-0.json
	-rm -f test_missinginput-0.json test_missinginput.json
