#include "analysis.h"
#include "sample.h"

#include "TChain.h"

void makeTables( analysis* myAnalysis );
void makeStack(  analysis* myAnalysis );
int ScanChain( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true);
int looperCR2lep( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true);
