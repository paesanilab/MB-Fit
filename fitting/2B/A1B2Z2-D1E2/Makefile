
CXX=g++
CXXFLAGS= -g -Wall -std=c++11 -O0 -m64 -I/opt/intel/mkl/include
LIBS = -lnetcdf -lgsl -lgslcblas -L/opt/intel/lib/intel64 -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

AR = /usr/bin/ar
OBJDIR = .
LIBDIR = ./
INCLUDE = -I./

FIT_OBJ = fit-utils.o coulomb.o electrostatics.o gammq.o io-xyz.o \
kvstring.o mon1.o mon2.o ps.o rwlsq.o wlsq.o stuff.o tang-toennies.o \
training_set.o ttm4.o poly_2b_A1B2Z2_D1E2_v1x.o \
x2b_A1B2Z2_D1E2_v1.o poly_2b_A1B2Z2_D1E2_v1.o \
dispersion.o poly_2b_A1B2Z2_D1E2.o 

EVAL_OBJ = fit-utils.o coulomb.o electrostatics.o gammq.o io-xyz.o \
kvstring.o mon1.o mon2.o ps.o rwlsq.o wlsq.o stuff.o tang-toennies.o \
training_set.o ttm4.o poly_2b_A1B2Z2_D1E2_v1x.o \
x2b_A1B2Z2_D1E2_v1x.o poly_2b_A1B2Z2_D1E2_v1.o \
dispersion.o poly_2b_A1B2Z2_D1E2.o

all: libfit.a libeval.a fit-2b eval-2b
# $(MAKE) libpot.a 2> err.txt > log.txt


libfit.a: $(addprefix $(OBJDIR)/, $(FIT_OBJ))
	$(AR) cru libfit.a $(addprefix $(OBJDIR)/, $(FIT_OBJ))

libeval.a: $(addprefix $(OBJDIR)/, $(EVAL_OBJ))
	$(AR) cru libeval.a $(addprefix $(OBJDIR)/, $(EVAL_OBJ))

fit-2b: fit-2b.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) $< $(LIBS) -lfit -o $@

eval-2b: eval-2b.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) $< $(LIBS) -leval -o $@

$(OBJDIR)/%.o: %.cpp $(OBJDIR)/.sentinel
	$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) -c $< $(LIBS) -o $@

$(OBJDIR)/.sentinel:
	mkdir -p $(OBJDIR)
	touch $@

clean:
	rm -rf $(addprefix $(OBJDIR)/, $(FIT_OBJ)) libfit.a libeval fit-2b eval-2b
