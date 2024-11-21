### MKL ###
# Information on how to link MKL can be found on 
# https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html

MKLROOT=/opt/intel/oneapi/mkl/2024.0/
MKL_CFLAGS=-DMKL_ILP64 -m64 -I${MKLROOT}/include
MKL_LDFLAGS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lm -ldl
### END MKL ###

### MPI ###
MPI_CFLAGS=$(shell mpicc --showme:compile)
MPI_LDFLAGS=$(shell mpicc --showme:link)
### END MPI ###

CC=gcc
CFLAGS=-m64 -Ofast -flto -march=native -funroll-loops $(MKL_CFLAGS) $(MPI_CFLAGS)
LDFLAGS=-lpthread -lm -ldl $(MKL_LDFLAGS) $(MPI_LDFLAGS)
EXEC=AGG

all: $(EXEC)

$(EXEC): main_mpi.o simul.o
		$(CC) -o $@ $^ $(LDFLAGS)

main_mpi.o: main_mpi.c simul.h
		$(CC) -o $@ -c $< $(CFLAGS)

simul.o: simul.c simul.h
		$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
		rm -f *.o

mrproper: clean
		rm -rf $(EXEC)
