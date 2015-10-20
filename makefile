# Define some commands

CXX         = g++  #simple. Just invokes the compiler.
CXXFLAGS    = -g -Wall -fPIC   #-g enables gdebug info. -Wall enables all warnings. -fPIC necessary for shared libraries.
ROOTCFLAGS  = $(shell root-config --cflags --libs)



runAll: ScanChain.o makeTables.o makeStack.o dataMCplotMaker.so sample.so analysis.so
	g++ $(CXXFLAGS) $(ROOTCFLAGS) runAll.cc ScanChain.o makeTables.o makeStack.o dataMCplotMaker.so sample.so analysis.so -Wl,-rpath,./ -o runAll


ScanChain.o: ScanChain.C
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -c ScanChain.C -o ScanChain.o

makeTables.o: makeTables.C
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -c makeTables.C -o makeTables.o

makeStack.o: makeStack.C
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -c makeStack.C -o makeStack.o

dataMCplotMaker.so: dataMCplotMaker.cc dataMCplotMaker.h PlotMakingTools.h
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -shared dataMCplotMaker.cc dataMCplotMaker.h PlotMakingTools.h -o dataMCplotMaker.so

sample.so: sample.cc
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -shared sample.cc -o sample.so

analysis.so: analysis.cc
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -shared analysis.cc -o analysis.so



.PHONY: clean

clean:
	rm -v *.o *.so *.d *.pcm
