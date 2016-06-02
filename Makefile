# Makefile

PROGRAM = main
SUB_OBJ = lu.o gen.o solve.o gauss.o rdft.o
OBJ     = main.o $(SUB_OBJ)
TESTER  = test.h blas_extern.h global.h
CC      = gcc
CFLAGS  = -Wall -O2 -std=c99
MOD_LIB = ./lapack_mod/mod_lu.a
HDR_DIR = -I/opt/OpenBLAS/include -I/opt/fftw/include
LIB_DIR = -L/opt/lapack -L/opt/OpenBLAS/lib -L/opt/fftw/lib
LIB     = -llapack -lopenblas -lfftw3 -lm -lgfortran
LIBLINK = $(HDR_DIR) $(LIB_DIR) $(LIB)
QM_CC   = gcc
QM_FLAG = -Wall -O2 -std=gnu99
QM_LIB  = -lquadmath


.PHONY: all
all: $(PROGRAM)

$(PROGRAM): $(OBJ) $(TESTER)
	$(CC) $(CFLAGS) -o main $(OBJ) $(MOD_LIB) $(LIBLINK) $(QM_LIB)

gen.o: gen.c $(TESTER)
	$(QM_CC) $(QM_FLAG) -c $< $(QM_LIB)

.c.o: $(TESTER)
	$(CC) $(CFLAGS) -c $< $(HDR_DIR)

.PHONY: clean
clean:
	rm $(PROGRAM) $(OBJ)
