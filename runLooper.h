#include "analysis.h"
#include "sample.h"

#include "TChain.h"

void makeTables( analysis* myAnalysis );
void makeStack(  analysis* myAnalysis );
void makeDataCards( analysis* myAnalysis );
void makeLostLepEstimate( analysis* srAnalysis, analysis* crAnalysis );
int ScanChain( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true);
int looperCR2lep( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true);
