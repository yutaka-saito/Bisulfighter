INCLUDE   = -I. -I/home/saito/local/include -I/usr/include -I/usr/local/include

LIBPATH   = -L. -L/home/saito/local/lib -L/usr/lib -L/usr/local/lib

LDLIBS    = -lm -lboost_program_options
#LDLIBS    = -lm -lboost_program_options -lboost_thread

DEFINE    = 
#DEFINE   += -DNDEBUG

CPPFLAGS += $(INCLUDE) $(DEFINE) -Wall -Wno-sign-compare -O3
#CPPFLAGS += -pthread
#CPPFLAGS += -static
#CPPFLAGS += -pg

LDFLAGS  += $(LIBPATH)

CXX ?= g++

PROGRAM1 = ComMet

OBJ1  = ComMet.o FrameworkComMet.o Utility.o \
	InputFormat.o Data.o \
	GlobalStatistics.o \
	ProbabilityModel.o NaiveModel.o DualModel.o \
	SlimProbabilityModel.o SlimNaiveModel.o SlimDualModel.o
HEAD1 = FrameworkComMet.hh LogSpace.hh Utility.hh \
	InputFormat.hh Data.hh \
	GlobalStatistics.hh \
	ProbabilityModel.hh NaiveModel.hh DualModel.hh \
	SlimProbabilityModel.hh SlimNaiveModel.hh SlimDualModel.hh

#------------- Make rules -----------------

.cc.o:
	$(CXX) $(CPPFLAGS) -c $< -o $@

.cpp.o:
	$(CXX) $(CPPFLAGS) -c $< -o $@

#------------- Commands -------------------

.PHONY: all clean

all: $(PROGRAM1) Makefile

clean:
	\rm -f *.o *~ $(PROGRAM1)
	chmod 644 *.cc *.hh

$(PROGRAM1): $(OBJ1)
	$(CXX) $(CPPFLAGS) -o $(PROGRAM1) $(OBJ1) $(LDFLAGS) $(LDLIBS)

$(OBJ1): $(HEAD1)
