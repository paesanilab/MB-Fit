CXX = icc
CXXFLAGS = -Wall -g -O0 -m64 -I/opt/intel/mkl/include 

LIBS = -lgsl -lgslcblas  -L/opt/intel/lib/intel64 -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

all: fit-1b-polynomials eval-1b-pot eval-1b-pot-wgrads 
test: eval-1b-pot-wgrads

fit-1b-polynomials: fit-1b-polynomials.cpp rwlsq.o wlsq.o x1b-v1.o training-set.o xyz-water-utils.o io-xyz.o stuff.o fit-utils.o 
	$(CXX) $(LIBS) $(CXXFLAGS) -o fit-1b-polynomials fit-1b-polynomials.cpp rwlsq.o wlsq.o x1b-v1.o training-set.o xyz-water-utils.o io-xyz.o stuff.o fit-utils.o	

eval-1b-pot: eval-1b-pot.cpp rwlsq.o wlsq.o x1b-v1.o training-set.o xyz-water-utils.o io-xyz.o stuff.o fit-utils.o
	$(CXX) $(LIBS) $(CXXFLAGS) -o eval-1b-pot eval-1b-pot.cpp rwlsq.o wlsq.o x1b-v1.o training-set.o xyz-water-utils.o io-xyz.o stuff.o fit-utils.o	

eval-1b-pot-wgrads: eval-1b-pot-wgrads.cpp  x1b-v1x.o poly-1b-v1x.o poly-1b-v1x-nogrd.o training-set.o xyz-water-utils.o io-xyz.o fit-utils.o
	$(CXX) $(LIBS) $(CXXFLAGS) -o eval-1b-pot-wgrads eval-1b-pot-wgrads.cpp poly-1b-v1x.o poly-1b-v1x-nogrd.o x1b-v1x.o training-set.o xyz-water-utils.o io-xyz.o fit-utils.o	

io-xyz.o: io-xyz.cpp io-xyz.h
rwlsq.o: rwlsq.cpp rwlsq.h
wlsq.o: wlsq.cpp wlsq.h
stuff.o: stuff.cpp stuff.h
training-set.o: training-set.cpp training-set.h
x1b-v1.o: x1b-v1.cpp x1b-v1.h stuff.o
x1b-v1x.o: x1b-v1x.cpp x1b-v1x.h poly-1b-v1x.cpp poly-1b-v1x-nogrd.cpp poly-1b-v1x.h
poly-1b-v1x.o: poly-1b-v1x.h poly-1b-v1x.cpp
poly-1b-v1x-nogrd.o: poly-1b-v1x-nogrd.cpp poly-1b-v1x.h	
xyz-water-utils.o: xyz-water-utils.cpp xyz-water-utils.h

fit-utils.o: fit-utils.h fit-utils.cpp
	$(CXX) -c $(CXXFLAGS) fit-utils.cpp   


clean:
	rm *.o fit-1b-polynomials eval-1b-pot eval-1b-pot-wgrads
