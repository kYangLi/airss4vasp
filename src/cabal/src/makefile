.DEFAULT:
	-touch $@
all: cabal
cabal.o: cabal.f90 niggli_reduce.o constants.o
	$(FC) $(FFLAGS) -c cabal.f90
niggli_reduce.o: ../../common/niggli_reduce.f90 constants.o
	$(FC) $(FFLAGS) -c ../../common/niggli_reduce.f90
constants.o: ../../common/constants.f90
constants.o: ../../common/constants.f90
	$(FC) $(FFLAGS) -c ../../common/constants.f90
SRC = cabal.f90 ../../common/niggli_reduce.f90 ../../common/constants.f90
OBJ = cabal.o niggli_reduce.o constants.o
clean: neat
	-rm -f cabal
neat:
	-rm -f $(OBJ) *.mod
cabal: $(OBJ) 
	$(LD) $(OBJ) -o cabal $(LDFLAGS)
