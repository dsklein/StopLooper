#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1.h"

#include "analysis.h"
#include "sample.h"

using namespace std;



void makeDataCards( analysis* myAnalysis ) {

  // Do the basic setup stuff

  vector<TString> bkgs = { "zNuNu", "dilep", "top1l", "W1l" }; // Eventually pull this from the analysis object
  const int nBkgs = bkgs.size();

  vector<TString> signals = myAnalysis->GetSignalLabels();
  const int nSigs = signals.size();
  if( nSigs > 1 ) {
	cout << "\nWarning in makeDataCards.cc! There are " << nSigs << " signals, when there should be just one!" << endl;
	cout << "  Using only the last signal (" << signals.at(signals.size()-1) << ").\n" << endl;
  }
  else if( nSigs < 1 ) {
	cout << "\nError in makeDataCards.cc: Need at least one signal sample!" << endl;
	return;
  }

  vector<TString> samples = bkgs;
  samples.insert( samples.begin(), signals.end()-1, signals.end() );
  const int nSamples = samples.size();

  vector<TString> sigRegions = myAnalysis->GetSigRegionsAll();
  const int nSigRegs = sigRegions.size();

  // Open files containing background yields and uncertainties
  TFile* uncertFile = new TFile( "uncertSR.root", "READ" );
  TFile* lostlepFile = new TFile( "bkgEstimates.root", "READ" );
  TH1D* h_lostLep = (TH1D*)lostlepFile->Get("lostLepBkg");
  
  //////////////////////////////////////////////////////////
  // Loop over signal regions, making a datacard for each...
  ///////////////////////////////////////////////////////////

  for( int bin=1; bin<=nSigRegs; bin++ ) {

	TH1D* h_yield = (TH1D*)uncertFile->Get( "evttype_"+sigRegions.at(bin-1) );

	// Do some acrobatics to send the output to a file...
	TString fileName = "datacards/datacard_"+sigRegions.at(bin-1)+".txt";
	cout << "Writing data card " << fileName << endl;

	FILE * outfile;
	outfile = fopen( fileName.Data(), "w" );


	////////////////////////////
	// Now print the data card

	fprintf( outfile,  "# Data card for signal region %d (%s)\n", bin, sigRegions.at(bin-1).Data() );
	fprintf( outfile,  "---\n" );
	fprintf( outfile,  "imax 1  number of channels\n" );

	fprintf( outfile,  "jmax %d  number of backgrounds (", nBkgs );
	for (TString bkgName : bkgs ) fprintf( outfile, "%s, ",  bkgName.Data() );   // Will usually be 4 (znunu, 2l, 1ltop, 1lw)
	fprintf( outfile,  ")\n" );

	fprintf( outfile,  "kmax %d  number of uncertainties\n", 2*nSamples-1 ); // For now, we've got statistical and systematic for each bkg, and stat for signal
	fprintf( outfile,  "---\n" );

	fprintf( outfile,  "# Now list the number of events observed (or zero if no data)\n" );
	fprintf( outfile,  "bin %d\n", bin );
	if( myAnalysis->HasData() ) fprintf( outfile,  "observation %d\n", int(h_yield->GetBinContent(1)) );
	else                        fprintf( outfile,  "observation 0\n" );
	fprintf( outfile,  "---\n" );

	fprintf( outfile,  "# Now list the expected events (i.e. Monte Carlo yield) for signal and all backgrounds in our particular bin\n" );
	fprintf( outfile,  "# Second 'process' line should be zero for signal, and a positive integer for each background\n" );

	fprintf( outfile,  "bin      ");
	for( TString sample : samples ) fprintf( outfile,  " %10i", bin );  // Print the signal region number once for each sample (signal and bkg)
	fprintf( outfile, "\n" );

	fprintf( outfile,  "process  ");
	for( TString sample : samples ) fprintf( outfile,  " %10s", sample.Data() );  // Print the name of each sample (sig & bkg)
	fprintf( outfile, "\n" );

	fprintf( outfile,  "process  ");
	for( int j=0; j<nSamples; j++ ) fprintf( outfile,  " %10i", j ); // Print a number for each sample (0=signal, positive integers for bkgs)
	fprintf( outfile, "\n" );



	// Print row with the yields for each process (signal, bkgs)
	fprintf( outfile,  "rate     ");
	double yield;

	for( int i=1; i<=nSamples; i++ ) {
	  if( i==3 ) yield = h_lostLep->GetBinContent(bin); // Pull 2l yield from lostLepton estimate histogram
	  else       yield = h_yield->GetBinContent(i+1);
	  fprintf( outfile,  " %10f", yield );
	}
	fprintf( outfile, "\n" );

	fprintf( outfile,  "---\n" );
	fprintf( outfile,  "# Now we list the independent sources of uncertainty (syst. and stat. error), and which samples they affect\n" );
	fprintf( outfile,  "---\n" );



	// Print out rows for statistical uncertainties
	for( int sampleIdx=0; sampleIdx<nSamples; sampleIdx++ ) {

	  char statname[25];
	  sprintf( statname, "%sStat%d", samples.at(sampleIdx).Data(), bin );	  
	  fprintf( outfile,   "%-18s  lnN ", statname );

	  double statErr = 1.0 + (h_yield->GetBinError(sampleIdx+2) / h_yield->GetBinContent(sampleIdx+2) );
	  if( sampleIdx==2 ) statErr = 1.0 + ( h_lostLep->GetBinError(bin) / h_lostLep->GetBinContent(bin) );  // Pull 2l stat error from lostLepton estimate histogram

	  for( int j=0; j<nSamples; j++ ) {
		if( j == sampleIdx )  fprintf( outfile, "  %8.6f  ", statErr);
		else fprintf( outfile,  "     -      " );
	  }
	  fprintf( outfile, "\n" );
	}


	// Print out rows for systematic uncertainties
	for( int sampleIdx = 1; sampleIdx<nSamples; sampleIdx++ ) { // For now, skip the signal sample (don't give it a systematic uncertainty)

	  char systname[25];
	  sprintf( systname, "%sSyst", samples.at(sampleIdx).Data() );	  
	  fprintf( outfile,   "%-18s  lnN ", systname );

	  double systErr = 1.3; // Flat 30% systematic for now
	  // if( sampleIdx = 4 ) systErr = 2.0; // 100% systematic on 1-lepton from top

	  for( int j=0; j<nSamples; j++ ) {
		if( j == sampleIdx )  fprintf( outfile, "  %8.6f  ", systErr);
		else fprintf( outfile,  "     -      " );
	  }
	  fprintf( outfile, "\n" );
	}


	fprintf( outfile,  "---\n" );
	fclose(outfile);

  } // End loop over signal regions

  uncertFile->Close();
  delete uncertFile;
}
