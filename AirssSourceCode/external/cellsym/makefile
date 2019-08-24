.DEFAULT:
	-touch $@

CELLSYM=cellsym-0.16a

INC=$(CELLSYM)/cellsym.h

SRC=$(CELLSYM)/cellsym.c $(CELLSYM)/cell_read.c $(CELLSYM)/periodic_table.c $(CELLSYM)/pdb_write.c $(CELLSYM)/cell_write.c $(CELLSYM)/pdb_read.c $(CELLSYM)/basis.c $(CELLSYM)/ident_sym.c $(CELLSYM)/shelx_read.c $(CELLSYM)/shelx_write.c

all: cellsym

cellsym:
	-wget -c --timeout=2 --tries=1 -nv http://www.tcm.phy.cam.ac.uk/sw/check2xsf/cellsym.tgz
	-tar -xzf cellsym.tgz
	-$(CC) $(CFLAGS) -L$(PREFIX)/lib -I$(PREFIX)/include -o $(CELLSYM)/cellsym $(INC) $(SRC) -lsymspg -lm

clean:
	-rm -fr $(CELLSYM)
