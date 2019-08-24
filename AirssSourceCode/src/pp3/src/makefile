.DEFAULT:
	-touch $@
all: pp3
pp3.o: pp3.f90 constants.o cell.o pp.o opt.o
	$(FC) $(FFLAGS) -c pp3.f90
cell.o: cell.f90 constants.o 
	$(FC) $(FFLAGS) -c cell.f90
pp.o: pp.f90 constants.o cell.o
	$(FC) $(FFLAGS) -c pp.f90
opt.o: opt.f90 constants.o cell.o pp.o
	$(FC) $(FFLAGS) -c opt.f90
constants.o: ../../common/constants.f90
	$(FC) $(FFLAGS) -c ../../common/constants.f90
SRC = opt.f90 pp.f90 cell.f90 pp3.f90 ../../common/constants.f90
OBJ = opt.o pp.o cell.o pp3.o constants.o
clean: neat
	-rm -f pp3
neat:
	-rm -f $(OBJ) *.mod
pp3: $(OBJ) 
	$(LD) $(OBJ) -o pp3
