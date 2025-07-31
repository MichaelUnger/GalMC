.PHONY: depend clean test

WITH_OPENMP = 1

LD := $(CXX)

DMWDIR = D2MW

ifndef HEALPIX
  HEALPIXDIR = /usr/include/healpix_cxx/
else
  HEALPIXDIR = $(HEALPIX)/src/cxx/$(HEALPIX_TARGET)
endif


CLING = rootcling
CXXFLAGS += -std=c++11 -O3  -Wl,--no-as-needed -Wno-unknown-pragmas -Wall -Wextra -fPIC -Wpedantic
CPPFLAGS += -I$(RUQIDIR)/include -I$(HEALPIXDIR) -I$(DMWDIR)/include
CPPFLAGS += -I$(shell root-config --incdir)
CPPFLAGS += -DGMCPATH=\"$(GMCROOT)\"

LDFLAGS  += -L$(RUQIDIR)/lib -lRUQI -lNE2001 -lYMW16 -lgfortran
LDFLAGS  += -Wl,--no-as-needed -fPIC $(shell root-config --ldflags)
LDFLAGS  += $(shell root-config --libs) -lEG
LDFLAGS  += -lboost_system -lboost_filesystem -lboost_program_options
LDFLAGS += -lgsl -lgslcblas
LDFLAGS  += -L$(HEALPIXDIR)/lib -lhealpix_cxx #-lcxxsupport -lc_utils -lcfitsio
SOFLAGS  += -fPIC -shared

ifeq ($(WITH_OPENMP), 1)
   CXXFLAGS += -D_WITH_OPENMP_ -fopenmp
   LDFLAGS  += -fopenmp
endif




LIBDIR := ./lib
SRCDIR := ./src
PRGDIR := ./prog
LIB := $(LIBDIR)/lib$(shell basename $(CURDIR)).so
OBJS := $(patsubst $(SRCDIR)/%.cc, $(SRCDIR)/%.o, $(wildcard $(SRCDIR)/*.cc))
OBJS += $(patsubst $(SRCDIR)/%LinkDef.h, $(SRCDIR)/%Dict.o, $(wildcard $(SRCDIR)/*LinkDef.h))
EXE = $(patsubst $(PRGDIR)/%.cxx, %, $(wildcard $(PRGDIR)/*.cxx))

all: depend $(LIB) $(EXE)

depend: Make-depend

$(SRCDIR)/%.o : $(SRCDIR)/%.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

$(SRCDIR)/%Dict.cc: $(shell grep -l ClassDef $(SRCDIR)/*.h) $(SRCDIR)/%LinkDef.h
	@(echo generating $@ dictionary)
	($(CLING) -f $@ -c $(CPPFLAGS) $^)
	@mkdir -p $(LIBDIR)
	@mv $*Dict_rdict.pcm $(LIBDIR)

$(LIB): $(OBJS)
	@mkdir -p $(LIBDIR)
	$(LD) $(SOFLAGS) $^ -o $@ $(LDFLAGS)

test: testGalMC
	./testGalMC --log_level=test_suite

%: $(PRGDIR)/%.cxx $(LIB)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -L$(LIBDIR) -lGalMC -lD2MW -lboost_unit_test_framework

Make-depend: $(wildcard *.h)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM -MG $^ > $@

clean:
	@rm -f $(LIB) $(OBJS) $(LIBDIR)/*.pcm Make-depend

ifneq ($(MAKECMDGOALS),clean)
-include Make-depend
endif
