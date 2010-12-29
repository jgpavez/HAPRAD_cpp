SHELL = /bin/bash

.DELETE_ON_ERROR:

.PHONY: all clean


ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
CERNLIBS    := -lpdflib804 -lmathlib -lphtools -lpacklib -lkernlib -lpawlib

CXX       := g++
CXXFLAGS  += -O2 -Wall -fPIC $(ROOTCFLAGS)
LD        = g++
LDFLAGS   = -O2 $(ROOTLDFLAGS) -lgfortran

INCLUDES  := -I..
LIBS      := -L../slib -lTRadCor $(ROOTLIBS) -lMathMore $(CERNLIBS)

##############################################################################

all: main_test

main_test: main_test.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

main_test.o: main_test.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -rf main_test main_test.o
