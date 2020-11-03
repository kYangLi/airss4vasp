## Please Change the parameters below...
# [tips] You may need to check: 
# [tips]  https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor 
# [tips]  to confirm the FFLAGS, FLIBS, CFLAGS, CLIBS, if you are using intel
# [tips]  compilers.
FC=ifort
FLIBS=${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
FFLAGS=-i8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include 

CC=icc
CLIBS=${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
CFLAGS=-DMKL_ILP64 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include


## DO NOT CHANGE BELOW ##
DFLAGS=-DCOMPAT

LD=$(FC)
LDFLAGS=$(FFLAGS) -L$(FLIBS)

LDC=$(CC)
LDCFLAGS=$(CFLAGS) -L$(CLIBS)

PREFIX=$(PWD)

export

all: internal external install

internal: cabal buildcell cryan genkp

external: spglib cellsym

cabal:
	(cd src/cabal/src; make)

buildcell:
	(cd src/buildcell/src; make)

cryan:
	(cd src/cryan/src; make)

genkp:
	(cd src/genkp; make)

spglib:
	(cd external/spglib; make)

cellsym:
	(cd external/cellsym; make)

install:
	(cp src/cabal/src/cabal bin/)
	(cp src/buildcell/src/buildcell bin/)
	(cp src/cryan/src/cryan bin/)
	(cp src/genkp/genkp bin/)
	(cp external/cellsym/cellsym-0.16a/cellsym bin/)
	@echo
	@echo 'Add '$(PWD)'/bin to your path by placing this line in your ~/.bashrc:'
	@echo 'export PATH='$(PWD)'/bin:$$PATH'
	@echo 'To update your path "source ~/.bashrc"'

neat_internal:
	(cd src/cabal/src; make clean)
	(cd src/buildcell/src; make clean)
	(cd src/cryan/src; make clean)
	(cd src/genkp; make clean)

neat_external:
	(cd external/cellsym; make clean)
	(cd external/spglib; make clean)

neat: neat_internal neat_external

clean: neat
	(rm -f bin/cabal bin/cryan bin/buildcell bin/cellsym)
