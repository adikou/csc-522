#------------------------------------------------------------------------------
# File: Makefile
#
# Note: This Makefile requires GNU make.
#
#------------------------------------------------------------------------------


CC = gcc
MPICC = mpicc
LIBS = -lm -lpthread
FLAGS = -O0
OMP_FLAGS = -fopenmp
SEQ_RB_SRC = seq-rb.c
MT_RB_SRC = mt-rb.c
DIST_RB_SRC = dist-rb.c
HYBRID_RB_SRC = hybrid-rb.c
BINARIES = seq-rb mt-rb dist-rb hybrid-rb

seq-rb : 
	$(CC) -o seq-rb $(FLAGS) $(SEQ_RB_SRC)

mt-rb :
	$(CC) -o mt-rb $(FLAGS) $(MT_RB_SRC) $(LIBS)

dist-rb : 
	$(MPICC) -o dist-rb $(FLAGS) $(DIST_RB_SRC)

hybrid-rb :
	$(MPICC) -o hybrid-rb $(OMP_FLAGS) $(FLAGS) $(HYBRID_RB_SRC)

.PHONY : clean

clean: 
	rm -f $(BINARIES) *.o *~ core *.dump

