#include "analysis.h"
#include "sample.h"

#include "TChain.h"

void makeTables( analysis* myAnalysis );
void makeStack(  analysis* myAnalysis );
int ScanChain( TChain* chain, std::string sampleName = "default", int nEvents = -1, bool fast = true);
