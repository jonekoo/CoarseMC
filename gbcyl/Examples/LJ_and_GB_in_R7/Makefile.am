
dist_data_DATA = input-0.json input-1.json

geometry.json : $(top_builddir)/src/create_crystal
	$(top_builddir)/src/create_crystal nx=40, ny=40, nz=12, boxtype=\"cylindrical\", radius=7.0, offset=1.0

ljgbr7-0.json ljgbr7-1.json : geometry.json input-0.json input-1.json $(top_srcdir)/src/update_json.py $(top_srcdir)/src/replace_with_points.py
	cp $(top_srcdir)/Examples/LJ_and_GB_in_R7/input-0.json $(top_srcdir)/Examples/LJ_and_GB_in_R7/input-1.json .
	mv input-0.json ljgbr7-0.json
	mv input-1.json ljgbr7-1.json
	$(PYTHON) $(top_srcdir)/src/update_json.py geometry.json ljgbr7-0.json ljgbr7-1.json
	$(PYTHON) $(top_srcdir)/src/replace_with_points.py ljgbr7-0.json ljgbr7-1.json

OUTPUT = restart-0.json restart-1.json output-0.json output-1.json

$(OUTPUT) : ljgbr7-0.json ljgbr7-1.json $(top_builddir)/src/ptgbcyl
	mpirun -np 2 $(top_builddir)/src/ptgbcyl -i ljgbr7-_I_.json 

mostlyclean-local : 
	-rm -f geometry.json
