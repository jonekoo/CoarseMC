
dist_data_DATA = input-0.json

geometry.json : $(top_builddir)/src/create_crystal
	@echo "\nCreating a hexagonal close-packed structure inside a \
	cylindrical cavity."
	$(top_builddir)/src/create_crystal nx=40, ny=40, nz=12, boxtype=\"cylindrical\", radius=7.0, offset=1.0


temperatures.json : $(top_srcdir)/src/geometric_progression.py
	@echo "\nCreating temperatures for the replicas."
	$(PYTHON) $(top_srcdir)/src/geometric_progression.py 1.98 2.0 2 > temperatures.json


ljgbr7-0.json ljgbr7-1.json : geometry.json input-0.json temperatures.json $(top_srcdir)/src/update_json.py $(top_srcdir)/src/replace_with_points.py $(top_srcdir)/src/scatter.py
	@echo "\nCreating input files from template."
	cp $(top_srcdir)/Examples/LJ_and_GB_in_R7/input-0.json ljgbr7-0.json
	cp $(top_srcdir)/Examples/LJ_and_GB_in_R7/input-0.json ljgbr7-1.json
	@echo "\nSetting temperatures to input files."
	$(PYTHON) $(top_srcdir)/src/scatter.py temperatures.json ljgbr7-_I_.json 2
	@echo "\nAdding coordinates to input files."
	$(PYTHON) $(top_srcdir)/src/update_json.py geometry.json ljgbr7-0.json ljgbr7-1.json
	@echo "\nReplacing randomly selected rods with points."
	$(PYTHON) $(top_srcdir)/src/replace_with_points.py 10 ljgbr7-0.json ljgbr7-1.json

example : $(OUTPUT)

OUTPUT = restart-0.json restart-1.json output-0.json output-1.json

$(OUTPUT) : ljgbr7-0.json ljgbr7-1.json $(top_builddir)/src/ptgbcyl
	@echo "\nCleaning previous output files out of the way."
	rm -f output-0.json output-1.json
	@echo "\nRunning a simulation with two replicas."
	@echo "This might take a while..."
	mpirun -np 2 $(top_builddir)/src/ptgbcyl -i ljgbr7-_I_.json 

mostlyclean-local : 
	-rm -f geometry.json
	-rm -f ljgbr7-0.json ljgbr7-1.json
	-rm -f output-0.json output-1.json
	-rm -f restart-0.json restart-1.json
	-rm -f temperatures.json
