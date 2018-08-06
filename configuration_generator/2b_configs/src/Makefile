CXX=icpc
MPICXX = mpicxx
CXXFLAGS= -g -Wall -O2 -std=c++11

AR = /usr/bin/ar
OBJDIR = ../obj
BINDIR = ../bin
LIBDIR = ../lib
LIBDIRS = -L../lib
INCLUDE =  -I/usr/include -I.

OBJ_TOOLS = random-rotation.o mt19937.o io-xyz.o io-utils.o read-xyz.o

LIBS = -ltools

EXES = libtools.a 2b_ts_rigid 2b_ts_flex

all: $(EXES) 
# $(MAKE) $(EXES) $(FITTING) 2> err.txt > log.txt

libtools.a: $(addprefix $(OBJDIR)/, $(OBJ_TOOLS))
	$(AR) cru libtools.a $(addprefix $(OBJDIR)/, $(OBJ_TOOLS))
	mkdir -p $(LIBDIR)
	mv libtools.a $(LIBDIR)

2b_ts_flex: 2b_ts_flex.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBDIRS) $< $(LIBS) -o $@
	mkdir -p $(BINDIR)
	mv 2b_ts_flex $(BINDIR)

2b_ts_rigid: 2b_ts_rigid.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBDIRS) $< $(LIBS) -o $@
	mkdir -p $(BINDIR)
	mv 2b_ts_rigid $(BINDIR)

$(OBJDIR)/%.o: %.cpp $(OBJDIR)/.sentinel
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

$(OBJDIR)/.sentinel:
	mkdir -p $(OBJDIR)
	touch $@

clean:
	rm -rf  2b_ts_flex 2b_ts_rigid $(OBJDIR) $(BINDIR) $(LIBDIR)

