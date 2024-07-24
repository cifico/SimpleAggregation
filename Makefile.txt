### MKL ###
# Information on how to link MKL can be found on 
# https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html

# Debian with package intel-mkl
#MKLLIBDIR=/usr/lib/x86_64-linux-gnu/
#MKL_CFLAGS=-DMKL_ILP64 -m64 -I/usr/include/mkl
#MKL_LDFLAGS=-Wl,--start-group $(MKLLIBDIR)/libmkl_intel_ilp64.a $(MKLLIBDIR)/libmkl_sequential.a $(MKLLIBDIR)/libmkl_core.a -Wl,--end-group -lm -ldl

# Custom installation in /home/lptl/intel/ (on gambit)
MKLROOT=/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/
# Or custom installation in /opt/intel-18/ (on ponyo or chihiro)
#MKLROOT=/opt/intel-18/compilers_and_libraries/linux/mkl/
MKL_CFLAGS=-DMKL_ILP64 -m64 -I${MKLROOT}/include
MKL_LDFLAGS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lm -ldl
### END MKL ###

### MPI ###
MPI_CFLAGS=$(shell mpicc --showme:compile)
MPI_LDFLAGS=$(shell mpicc --showme:link)
### END MPI ###

CC=gcc
CFLAGS=-m64 -Ofast -flto -march=native -funroll-loops -DLARGE_NOBS $(MKL_CFLAGS) $(MPI_CFLAGS)
#CC=icc
#CFLAGS=-W -Wall -xSSE2 -axAVX,CORE-AVX2,CORE-AVX512 -Ofast $(MKL_CFLAGS) $(MPI_CFLAGS)
LDFLAGS=-lpthread -lm -ldl $(MKL_LDFLAGS) $(MPI_LDFLAGS)
EXEC=SEPfsaveGC

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
