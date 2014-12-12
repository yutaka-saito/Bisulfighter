INCLUDE   = -I. -I/usr/include -I/usr/local/include 

LIBPATH   = -L. -L/usr/lib -L/usr/local/lib 

LDLIBS    = -lm -lboost_program_options
#LDLIBS    = -lm -lboost_program_options -lboost_thread

DEFINE    = 
DEFINE   += -DNDEBUG

CPPFLAGS += $(INCLUDE) $(DEFINE) -Wall -Wno-sign-compare -O3
#CPPFLAGS += -pthread 
#CPPFLAGS += -pg 

LDFLAGS  += $(LIBPATH)

CXX ?= g++

PROGRAM1 = ComMet

OBJ1  = ComMet.o FrameworkComMet.o Utility.o \
	InputFormat.o Data.o \
	ProbabilityModel.o NaiveModel.o CGIModel.o \
	SlimProbabilityModel.o SlimNaiveModel.o SlimCGIModel.o
HEAD1 = FrameworkComMet.h LogSpace.hpp Utility.h \
	InputFormat.h Data.h \
	ProbabilityModel.h NaiveModel.h CGIModel.h \
	SlimProbabilityModel.h SlimNaiveModel.h SlimCGIModel.h

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
	chmod 644 *.cc *.h *.hpp

$(PROGRAM1): $(OBJ1)
	$(CXX) $(CPPFLAGS) -o $(PROGRAM1) $(OBJ1) $(LDFLAGS) $(LDLIBS)

$(OBJ1): $(HEAD1)
