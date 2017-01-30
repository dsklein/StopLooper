#include "analysis.h"
#include "sample.h"

#include "TChain.h"

void makeTables( analysis* myAnalysis );
void makeStack(  analysis* myAnalysis );
void makeDataCards( analysis* srAnalysis, analysis* sigAnalysis = NULL, analysis* lostlepAnalysis = NULL, analysis* onelepwAnalysis = NULL );
void makeLostLepEstimate( analysis* srAnalysis, analysis* crAnalysis );
void make1lWEstimate( analysis* srAnalysis, analysis* crAnalysis );
void makeSignalEstimate( analysis* sigAnalysis, analysis* genmetAnalysis = NULL );
int ScanChain( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true);
int looperCR2lep( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true);
int looperCR0b( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true);
