F90=gfortran
FLAGS=-std=f95 -pedantic -Wall


OBJS = nrtype.o particle.o orientational_ordering.o state_reader.o utils.o orientation_parameter_wrapper.o cylindrical_orientation_wrapper.o cylindrical_orientation.o translational_ordering.o tau1_wrapper.o 
MODS = nrtype.mod particle.mod orientational_ordering.mod utils.mod state_reader.mod orientation_parameter_wrapper.mod cylindrical_orientation_wrapper.mod cylindrical_orientation.mod translational_ordering.mod tau1_wrapper.mod 


all : orientation_parameter cylindrical_orientation_tensor tau1_parameter

orientation_parameter : $(OBJS) orientation_parameter.o
	$(F90) $(OBJS) orientation_parameter.o $(FLAGS) -o orientation_parameter

cylindrical_orientation_tensor : $(OBJS) cylindrical_orientation_tensor.o
	$(F90) $(OBJS) cylindrical_orientation_tensor.o $(FLAGS) -o cylindrical_orientation_tensor

tau1_parameter : $(OBJS) tau1_parameter.o
	$(F90) $(OBJS) tau1_parameter.o $(FLAGS) -o tau1_parameter

%.o %.mod : %.f90
	$(F90) -c $(FLAGS) $<


tau1_wrapper.o tau1_wrapper.mod : translational_ordering.mod orientational_ordering.mod nrtype.mod particle.mod 
cylindrical_orientation_wrapper.mod cylindrical_orientation_wrapper.o : particle.mod nrtype.mod cylindrical_orientation.mod
cylindrical_orientation.mod cylindrical_orientation.o : particle.mod nrtype.mod orientational_ordering.mod
nrtype.o nrtype.mod : nrtype.f90
particle.o particle.mod : particle.f90 nrtype.mod utils.mod
utils.o utils.mod : utils.f90 nrtype.mod
orientational_ordering.o orientational_ordering.mod : orientational_ordering.f90 particle.mod 
orientation_parameter_wrapper.mod orientation_parameter_wrapper.o : orientational_ordering.mod particle.mod nrtype.mod
state_reader.o state_reader.mod : state_reader.f90 particle.mod nrtype.mod
orientation_parameter.o : orientation_parameter_wrapper.mod state_reader.mod particle.mod nrtype.mod
cylindrical_orientation_tensor.o : cylindrical_orientation_wrapper.mod state_reader.mod particle.mod nrtype.mod
tau1_parameter.o : tau1_wrapper.mod state_reader.mod particle.mod nrtype.mod


clean :
	rm cylindrical_orientation_tensor.o orientation_parameter.o tau1_parameter.o $(OBJS) $(MODS)

