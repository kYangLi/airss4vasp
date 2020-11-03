FC=gfortran   # GCC family only
FFLAGS=-O0 -g # Do not optimise
#FFLAGS=-O0 -g -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -fall-intrinsics -ffpe-trap=invalid,zero,overflow -fbounds-check -Wuninitialized

CC=gcc # GCC family only
CFLAGS=-lm

LD=$(FC)
LDFLAGS=-llapack


GCC_VER_GTE46 := $(shell echo `$(FC) -dumpfullversion -dumpversion | \
                   cut -f1-2 -d.` \>= 4.6 | bc )
ifeq ($(GCC_VER_GTE46),0)
DFLAGS=-DCOMPAT
endif

LDC=$(CC)
LDCFLAGS=$(CFLAGS)

PREFIX=$(PWD)

export

all: internal external

internal: cabal buildcell cryan

external: spglib cellsym

cabal:
	(cd src/cabal/src; make)

buildcell:
	(cd src/buildcell/src; make)

cryan:
	(cd src/cryan/src; make)

spglib:
	(cd external/spglib; make)

cellsym:
	(cd external/cellsym; make)

install:
	(cp src/cabal/src/cabal bin/)
	(cp src/buildcell/src/buildcell bin/)
	(cp src/cryan/src/cryan bin/)
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
	(rm -f bin/cabal bin/cryan bin/buildcell bin/cellsym)
