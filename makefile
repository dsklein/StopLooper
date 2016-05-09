# Define some commands

CXX         = g++  #simple. Just invokes the compiler.
CXXFLAGS    = -g -Wall -fPIC   #-g enables gdebug info. -Wall enables all warnings. -fPIC necessary for shared libraries.
ROOTCFLAGS  = $(shell root-config --cflags --libs)



runLooper: runLooper.cc runLooper.h ScanChain.o makeTables.o makeStack.o libdataMCplotMaker.so libsample.so libanalysis.so libsfManager.so
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -lGenVector runLooper.cc ScanChain.o makeTables.o makeStack.o -L. -Wl,-rpath,./ -ldataMCplotMaker -lsample -lanalysis -lsfManager -o runLooper
#GenVector seems to be necessary to take a LorentzVector invariant mass. The error was "undefined reference to Math::GenVector::Throw()"


ScanChain.o: ScanChain.C
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -c ScanChain.C -o ScanChain.o

makeTables.o: makeTables.cc
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -c makeTables.cc -o makeTables.o

makeStack.o: makeStack.cc
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -c makeStack.cc -o makeStack.o

libdataMCplotMaker.so: dataMCplotMaker.cc dataMCplotMaker.h PlotMakingTools.h
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -shared dataMCplotMaker.cc dataMCplotMaker.h PlotMakingTools.h -o libdataMCplotMaker.so

libsample.so: sample.cc sample.h
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -shared sample.cc -o libsample.so

libanalysis.so: analysis.cc analysis.h
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -shared analysis.cc -o libanalysis.so

libsfManager.so: sfManager.cc sfManager.h btagsf/BTagCalibrationStandalone.cc
	g++ $(CXXFLAGS) $(ROOTCFLAGS) -shared sfManager.cc btagsf/BTagCalibrationStandalone.cc -o libsfManager.so



.PHONY: clean

clean:
	rm -v -f *.o *.so *.d *.pcm runLooper
