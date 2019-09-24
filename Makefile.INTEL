FC=ifort
FFLAGS=-O0 -g # Do not optimise
#FFLAGS= -O0 -g -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -fall-intrinsics -ffpe-trap=invalid,zero,overflow -fbounds-check -Wuninitialized

CC=icc
CFLAGS=

# [tips] You may need to check: 
# [tips] https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor 
# [tips] to confirm the LDFLAGS, if you are using intel compilers.
LDFLAGS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
LD=$(FC)

LDC=$(CC)
LDCFLAGS=-lm

# 'DFLAGS=-DCOMPAT' is necessary when it comes the error:
#   [COMPILER_VERSION] and [COMPILER_OPTIONS] are invalid in cell.f90.
#DFLAGS=-DCOMPAT

PREFIX=$(PWD)

export

all: internal external

internal: cabal buildcell cryan kpgen 

external: spglib cellsym

cabal:
	(cd src/cabal/src; make)

buildcell:
	(cd src/buildcell/src; make)

cryan:
	(cd src/cryan/src; make)

kpgen:
	(cd src/kpgen; make)

spglib:
	(cd external/spglib; make)

cellsym:
	(cd external/cellsym; make)

install:
	(cp src/cabal/src/cabal bin/)
	(cp src/buildcell/src/buildcell bin/)
	(cp src/cryan/src/cryan bin/)
	(cp src/kpgen/kpgen bin/)
	(cp external/cellsym/cellsym-0.16a/cellsym bin/)
	@echo
	@echo 'Add '$(PWD)'/bin to your path by placing this line in your ~/.bashrc:'
	@echo 'export PATH='$(PWD)'/bin:$$PATH'
	@echo 'To update your path "source ~/.bashrc"'

neat_internal:
	(cd src/cabal/src; make clean)
	(cd src/buildcell/src; make clean)
	(cd src/cryan/src; make clean)

neat_external:
	(cd external/cellsym; make clean)
	(cd external/spglib; make clean)

neat: neat_internal neat_external

clean: neat
	(rm -f bin/pp3 bin/cabal bin/cryan bin/buildcell bin/symmol bin/cellsym)
