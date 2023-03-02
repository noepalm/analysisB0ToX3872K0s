ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA -lRooFit -lRooFitCore -lMinuit

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g -lGenVector # con errori tipo 'undefined reference to `ROOT::Math::GenVector::Throw(char const*)'
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew , $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)B0toX3872K0s_base.o: $(INCLUDEDIR)src/B0toX3872K0s_base.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)B0toX3872K0s_base.o $<
$(OUTLIB)HLTapply.o: $(INCLUDEDIR)src/HLTapply.C 
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HLTapply.o $<

# ==================== HLT EMULATION =============================================
HLTemul:$(INCLUDEDIR)src/X3872Application.cc\
		$(OUTLIB)B0toX3872K0s_base.o\
		$(OUTLIB)HLTapply.o
		$(CXX) $(CXXFLAGS) -ldl -o X3872Application $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
HLTemul.clean:
		rm -f HLTemul 

hlt: HLTemul

clean:
		rm -f $(OUTLIB)*.o
		rm -f HLTemul


