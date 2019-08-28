OS            := $(shell uname)

#ifeq ($(OS), Linux)
CXX           := g++
DEBUG         := -g
#endif

CXXFLAGS      := $(DEBUG) -ansi -Wall -Wextra -m64 -O3 -std=c++11

SOURCEDIR     := ../src
SOURCES       := $(filter-out $(SOURCEDIR)/B0KstMuMu.cc,$(wildcard $(SOURCEDIR)/*.cc)) $(SOURCEDIR)/RooDataHist.cxx

INCLUDEDIR    := ../interface
INCLUDES      := -I $(INCLUDEDIR)

XERCESINFLAGS := -I $(HOME)/xerces/include
BOOSTINFLAGS  := -I $(HOME)/boost
GSLINCFLAGS   := -I $(HOME)/gsl/include

XERCESFLAGS   := -L $(HOME)/xerces/lib -lxerces-c
BOOSTFLFLAGS  := -L $(HOME)/boost/stage/lib
GSLFLAGS      := -L $(HOME)/gsl/lib -lgsl

ROOTFLAGS     := $(shell root-config --glibs)
ROOFITFLAGS   := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore -lMathMore -lMinuit -lFoam
#compile class
ROOTCINT := $(shell which rootcint)

CLASS2       := RooBernsteinEffi
CLASSDICT2   := $(CLASS2)Dictionary.cxx
LIBS2 := $(CLASS2).cxx  $(CLASSDICT2)
$(CLASSDICT2): $(CLASS2).h $(CLASS2)LinkDef.h
	$(ROOTCINT) -f $@ -c $^

CLASSDICT3  := AngDict.cc
$(CLASSDICT3):  $(INCLUDEDIR)/PdfRT.h  $(INCLUDEDIR)/PdfWT.h $(INCLUDEDIR)/PdfSigAng.h $(INCLUDEDIR)/DecayRate.h 
	$(ROOTCINT) -f $@ -c $^


ExtractYield : $(SOURCEDIR) $(INCLUDEDIR) ExtractYield.cc
	$(CXX) $(CXXFLAGS) -o $@  $(ROOFITFLAGS) $(INCLUDES) $(SOURCES) $(CLASSDICT3)  -D"ROOFIT" $@.cc
Moment : Moment.cc
	$(CXX) $(CXXFLAGS) -o $@  $(ROOTFLAGS) $(ROOFITFLAGS) $(LIBS2)  -D"ROOFIT" $@.cc
GenerateCocktail: $(SOURCEDIR) $(INCLUDEDIR) GenerateCocktail.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOFITFLAGS) $(INCLUDES) $(SOURCES)   -D"ROOFIT" $@.cc
MoMToy : MoMToy.cc
	$(CXX) $(CXXFLAGS) -o $@  $(ROOTFLAGS) $(ROOFITFLAGS)  -D"ROOFIT" $@.cc



LOCALCC := $(wildcard *.o)

#cleaning options
.PHONY: clean
clean:
	$(RM) $(LOCALCC:%.cc=%)  $(CLASSDICT1) $(CLASSDICT2)  Macros_C.*
