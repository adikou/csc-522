#------------------------------------------------------------------------------
# File: Makefile
#
# Note: This Makefile requires GNU make.
#
#------------------------------------------------------------------------------

MPICC = mpicc
PROFILER_SRC = profiler.c
BINARIES = libprofiler.a

profiler : 
	$(MPICC) -c $(PROFILER_SRC) -o $(BINARIES)

.PHONY : clean

clean: 
	rm -f $(BINARIES) *.o *~ core *.dump

