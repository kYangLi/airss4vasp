.DEFAULT:
	-touch $@
all: buildcell
opt.o: opt.f90 constants.o cell.o pp.o
	$(FC) $(FFLAGS) -c opt.f90
pp.o: pp.f90 constants.o cell.o symmetry.o
	$(FC) $(FFLAGS) -c pp.f90
build.o: build.f90 constants.o cell.o symmetry.o rng.o
	$(FC) $(FFLAGS) -c build.f90
buildcell.o: buildcell.f90 constants.o cell.o symmetry.o rng.o build.o pp.o opt.o
	$(FC) $(FFLAGS) -c buildcell.f90
cell.o: cell.f90 constants.o symmetry.o rng.o
	$(FC) $(FFLAGS) $(DFLAGS) -cpp -c cell.f90
constants.o: ../../common/constants.f90
	$(FC) $(FFLAGS) -c ../../common/constants.f90
rng.o: ../../common/rng.f90 constants.o
	$(FC) $(FFLAGS) -c ../../common/rng.f90
symmetry.o: symmetry.f90 constants.o
	$(FC) -O0 -c symmetry.f90
niggli_reduce.o: ../../common/niggli_reduce.f90 constants.o
	$(FC) $(FFLAGS) -c ../../common/niggli_reduce.f90
SRC = ../../common/niggli_reduce.f90 symmetry.f90 ../../common/rng.f90 cell.f90 buildcell.f90 build.f90 ../../common/constants.f90 pp.f90 opt.f90
OBJ = niggli_reduce.o symmetry.o rng.o cell.o buildcell.o build.o constants.o pp.o opt.o
clean: neat
	-rm -f buildcell
neat:
	-rm -f $(OBJ) *.mod
buildcell: $(OBJ) 
	$(LD) $(OBJ) -o buildcell
