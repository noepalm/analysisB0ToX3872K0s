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
$(OUTLIB)MVAoptimizer.o: $(INCLUDEDIR)src/MVAoptimizer.C $(INCLUDEDIR)include/LumiConstants.h
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MVAoptimizer.o $<
$(OUTLIB)RootPlots.o: $(INCLUDEDIR)tools/RootPlots.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RootPlots.o $<
# ==================== HLT EMULATION =============================================
HLTemul:$(INCLUDEDIR)src/X3872Application.cc\
		$(OUTLIB)B0toX3872K0s_base.o\
		$(OUTLIB)HLTapply.o
		$(CXX) $(CXXFLAGS) -ldl -o X3872Application $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
HLTemul.clean:
		rm -f HLTemul 

MVAopt:$(INCLUDEDIR)src/MVAoptimization.cc\
	   $(OUTLIB)MVAoptimizer.o
	   $(CXX) $(CXXFLAGS) -ldl -o MVAoptimization $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
MVAopt.clean:
	   rm -f MVAopt	

MVAcut:$(INCLUDEDIR)src/MVAcutter.cc\
	   $(OUTLIB)MVAoptimizer.o\
	   $(OUTLIB)RootPlots.o
	   $(CXX) $(CXXFLAGS) -ldl -o MVAcutter $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
MVAcut.clean:
		rm -f MVAcut

hlt: HLTemul
opt: MVAopt
cut: MVAcut

clean:
		rm -f $(OUTLIB)*.o
		rm -f HLTemul
		rm -f MVAopt
		rm -f MVAcut


