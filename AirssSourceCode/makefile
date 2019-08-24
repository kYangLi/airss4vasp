FC=gfortran   # GCC family only
FFLAGS=-O0 -g # Do not optimise

#FFLAGS= -O0 -g -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -fall-intrinsics -ffpe-trap=invalid,zero,overflow -fbounds-check -Wuninitialized

CC=gcc # GCC family only
CFLAGS=

LDFLAGS=-llapack
LD=$(FC)
PREFIX=$(PWD)

GCC_VER_GTE46 := $(shell echo `$(FC) -dumpfullversion -dumpversion | cut -f1-2 -d.` \>= 4.6 | bc )
ifeq ($(GCC_VER_GTE46),0)
DFLAGS=-DCOMPAT
endif

export

all: internal external

internal: pp3 cabal buildcell cryan 

external: symmol spglib cellsym

pp3:
	(cd src/pp3/src; make)

cabal:
	(cd src/cabal/src; make)

buildcell:
	(cd src/buildcell/src; make)

cryan:
	(cd src/cryan/src; make)

spglib:
	(cd external/spglib; make)

symmol:
	(cd external/symmol; make)

cellsym:
	(cd external/cellsym; make)

install:
	(cp src/pp3/src/pp3 bin/)
	(cp src/cabal/src/cabal bin/)
	(cp src/buildcell/src/buildcell bin/)
	(cp src/cryan/src/cryan bin/)
	(cp external/symmol/symmol bin/)
	(cp external/cellsym/cellsym-0.16a/cellsym bin/)
	@echo
	@echo 'Add '$(PWD)'/bin to your path by placing this line in your ~/.bashrc:'
	@echo 'export PATH='$(PWD)'/bin:$$PATH'
	@echo 'To update your path "source ~/.bashrc"'
	@([ -d .hg ] && (printf "echo ";head -1 VERSION | tr -d '\n';printf " build ";hg log -r . --template "{node|short}+ {branches} {date|rfc822date}") > bin/airss_version) || echo " "
	@chmod u+x bin/airss_version

check:
	(bash bin/check_airss)

neat_internal:
	(cd src/pp3/src; make clean)
	(cd src/cabal/src; make clean)
	(cd src/buildcell/src; make clean)
	(cd src/cryan/src; make clean)

neat_external:
	(cd external/symmol; make clean)
	(cd external/cellsym; make clean)
	(cd external/spglib; make clean)

neat: neat_internal neat_external

clean: neat
	(rm -f bin/pp3 bin/cabal bin/cryan bin/buildcell bin/symmol bin/cellsym)

dist: clean
	rm -fr .check external/cellsym/cellsym*.tgz external/spglib/spglib-*.tar.gz external/symmol/symmol.zip
	rm -fr include/* lib/*
	find examples/ -name "*-*.*" | xargs rm	-f
	find examples/ -name ".spawnpid*" | xargs rm -f
	#tar -czf ../airss-`date "+%d%m%Y"`.tgz  --exclude=".*" -C .. airss
	@echo 'For release use "hg archive -t tgz ../airss-x.y.z.tgz"'
