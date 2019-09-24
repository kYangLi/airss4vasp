.DEFAULT:
	-touch $@
all: cryan
cryan.o: cryan.f90 rng.o algor.o gp.o modularity.o constants.o
	$(FC) $(FFLAGS) -c cryan.f90
modularity.o: modularity.f90 rng.o constants.o
	$(FC) $(FFLAGS) -c modularity.f90
gp.o: gp.f90 algor.o constants.o
	$(FC) $(FFLAGS) -c gp.f90
constants.o: ../../common/constants.f90
	$(FC) $(FFLAGS) -c ../../common/constants.f90
algor.o: algor.f90 constants.o
	$(FC) $(FFLAGS) -c algor.f90
rng.o: ../../common/rng.f90 constants.o
	$(FC) $(FFLAGS) -c ../../common/rng.f90
SRC = cryan.f90 modularity.f90 gp.f90 algor.f90 ../../common/rng.f90 ../../common/constants.f90
OBJ = cryan.o modularity.o gp.o algor.o rng.o constants.o
clean: neat
	-rm -f cryan
neat:
	-rm -f $(OBJ) *.mod
cryan: $(OBJ) 
	$(LD) $(OBJ) -o cryan $(LDFLAGS)
