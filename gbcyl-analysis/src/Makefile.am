bin_PROGRAMS = orientation_parameter cylindrical_orientation tau1

orientation_parameter_SOURCES = nrtype.f90 utils.f90 particle.f90 \
state_reader.f90 orientational_ordering.f90 orientation_parameter_wrapper.f90 \
orientation_parameter.f90

cylindrical_orientation_SOURCES = nrtype.f90 utils.f90 particle.f90 \
state_reader.f90 orientational_ordering.f90 \
cylindrical_orientation.f90 cylindrical_orientation_wrapper.f90 \
cylindrical_orientation_tensor.f90

tau1_SOURCES = nrtype.f90 nrutil.f90 mnbrak.f90 nr.f90 tau1_module.f90 

MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)

