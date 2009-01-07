OBJS=main.o mpicelllist.o linkedlist.o particle.o nrtype.o io.o particlearray.o
MODS=mpicelllist.mod linkedlist.mod particle.mod nrtype.mod io.mod mpi.mod mpi_sizeofs.mod particlearray.mod

F90=mpif90

pargbcyl : $(OBJS) $(MODS)
	$(F90) -o pargbcyl $(OBJS)

%.o %.mod : %.f90
	$(F90) -c $<

main.o : mcengine.mod mpicelllist.mod particle.mod nrtype.mod io.mod particlearray.mod
mpicelllist.o mpicelllist.mod : linkedlist.mod particle.mod nrtype.mod mpi.mod mpi_sizeofs.mod particlearray.mod
particlearray.o particlearray.mod : particle.mod 
particle.o particle.mod : nrtype.mod 
linkedlist.o linkedlist.mod : particle.mod particlearray.mod
nrtype.o nrtype.mod : 
mpi.mod :
mpi_sizeofs.mod :
io.o io.mod : nrtype.mod 

clean :
	rm $(OBJS) $(MODS)
 
