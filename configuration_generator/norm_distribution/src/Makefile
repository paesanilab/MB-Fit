CXX = g++ #icpc
FC = gfortran #ifort
CXXFLAGS = -O2 -std=c++11
#FFLAGS = -O2
FFLAGS = -O0 -Wall -Wextra -pedantic -fbounds-check #-O0 -check bounds -g -heap-arrays -traceback

MKLFLAGS = #-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

LD = $(CXX)
LDFLAGS	= $(CXXFLAGS)


all: generate_configs_normdistrbn

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^


generate_configs_normdistrbn: generate_configs_normdistrbn.o sobol.o sobol_stdnormal.o random_stdnormal.o 
	$(FC) -o $@ $^ $(FFLAGS) $(MKLFLAGS) -lstdc++

clean:
	rm -f *.o *.mod generate_configs_normdistrbn

