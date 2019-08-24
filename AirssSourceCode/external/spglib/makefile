.DEFAULT:

all: get libsymspg.a

SPGLIB=spglib-1.9.4

OBJ = $(SPGLIB)/src/arithmetic.o $(SPGLIB)/src/cell.o $(SPGLIB)/src/debug.o $(SPGLIB)/src/delaunay.o $(SPGLIB)/src/hall_symbol.o $(SPGLIB)/src/kgrid.o \
    $(SPGLIB)/src/kpoint.o $(SPGLIB)/src/mathfunc.o $(SPGLIB)/src/niggli.o $(SPGLIB)/src/pointgroup.o $(SPGLIB)/src/primitive.o \
    $(SPGLIB)/src/refinement.o $(SPGLIB)/src/site_symmetry.o $(SPGLIB)/src/sitesym_database.o $(SPGLIB)/src/spacegroup.o \
    $(SPGLIB)/src/spg_database.o $(SPGLIB)/src/spglib.o $(SPGLIB)/src/spin.o $(SPGLIB)/src/symmetry.o

get:
	-wget -c --timeout=2 --tries=1 -nv http://www.check2xsf.org.uk/$(SPGLIB).tar.gz
	-tar -xzf $(SPGLIB).tar.gz

libsymspg.a: $(OBJ)
	-ar -cr $(PREFIX)/lib/libsymspg.a $(OBJ)
	-mkdir -p $(PREFIX)/include/spglib
	-cp $(SPGLIB)/src/spglib.h $(PREFIX)/include/spglib/spglib.h

clean:
	-rm -fr $(SPGLIB)
