.DEFAULT:
	-touch $@
all: genkp
genkp.o: genkp.c
	$(CC) $(CFLAGS) -c genkp.c
SRC = genkp.c
OBJ = genkp.o
clean: neat
	-rm -f genkp
neat:
	-rm -f $(OBJ)
genkp: $(OBJ) 
	$(LDC) $(OBJ) -o genkp $(LDCFLAGS)