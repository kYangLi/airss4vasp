.DEFAULT:
	-touch $@
all: kpgen
kpgen.o: kpgen.c
	$(CC) $(CFLAGS) -c kpgen.c
SRC = kpgen.c
OBJ = kpgen.o
kpgen: $(OBJ) 
	$(LDC) $(OBJ) -o kpgen $(LDCFLAGS)
