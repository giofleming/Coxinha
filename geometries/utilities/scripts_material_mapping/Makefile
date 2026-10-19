ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
                                                                                           
CXX = g++ -std=c++20
CXXFLAGS = -O4 $(ROOTCFLAGS)
DEPEND = $(CXX) -MM
DEL = rm -rf
                                                                                           
LD = $(CXX)
LDFLAGS = $(ROOTLIBS)
# LDFS = -lstdc++fs
                                                                                           
SOURCES           := $(wildcard *.cc)
OBJECTS           := $(SOURCES:.cc=.o)
EXESRC            := $(wildcard *.cpp)
EXECUTABLES       := $(EXESRC:.cpp=.exe)
all: $(EXECUTABLES)
                                                                                           
%.exe : %.cpp $(OBJECTS); ${LD} $(CXXFLAGS) $< $(OBJECTS) $(LDFLAGS) $(LDFS) -o $@ 
                                                                                           
%.o : %.cc; $(CXX) $(CXXFLAGS) -c $< -o $@
                                                                                           
%.d : %.cc; $(DEPEND) $(CXXFLAGS) $< > $@
                                                                                           
clean: ; $(DEL) *.d; $(DEL) *.o; $(DEL) *.exe
                                                                                           
clear: ; $(DEL) *~
                                                                                           
include $(SOURCES:.cc=.d)
