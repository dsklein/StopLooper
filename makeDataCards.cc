#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1.h"

#include "analysis.h"
#include "sample.h"

using namespace std;



void makeDataCards( analysis* myAnalysis ) {

  // Do the basic setup stuff

  vector<TString> bkgs = myAnalysis->GetBkgLabels();
  const int nBkgs = bkgs.size();

  vector<TString> signals = myAnalysis->GetSignalLabels();
  const int nSigs = signals.size();
  if ( nSigs != 1 ) {
	cout << "\nWarning in makeDataCards.cc! There are " << nSigs << " signals, when there should be just one!" << endl;
	cout << "  Using only the last signal (" << signals.at(signals.size()-1) << ").\n" << endl;
	// return;
  }

  vector<TString> samples = bkgs;
  samples.insert( samples.begin(), signals.end()-1, signals.end() );
  const int nSamples = samples.size();

  vector<TString> sigRegions = myAnalysis->GetSigRegionsAll();
  const int nSigRegs = sigRegions.size();

  
  //////////////////////////////////////////////////////////
  // Loop over signal regions, making a datacard for each...
  ///////////////////////////////////////////////////////////

  for( int bin=1; bin<=nSigRegs; bin++ ) {

	// Do some fancy trickery to send the output to a file...
	TString fileName = "datacards/datacard_"+sigRegions.at(bin-1)+".txt";
	cout << "Writing data card " << fileName << endl;

	FILE * outfile;
	outfile = fopen( fileName.Data(), "w" );

	// freopen( fileName.Data(), "w", stdout );

	// ofstream outfile;
	// streambuf* restore;
	// outfile.open( Form("combine_files/datacard_%s.txt", sigRegions.at(bin-1).Data()) );
	// restore = cout.rdbuf();
	// cout.rdbuf(outfile.rdbuf());

	// Open file containing uncertainty histograms
	TFile* histfile = new TFile( "uncertainties.root", "READ" );

	////////////////////////////
	// Now print the data card

	fprintf( outfile,  "# Data card for signal region %d (%s)\n", bin, sigRegions.at(bin-1).Data() );
	fprintf( outfile,  "---\n" );
	fprintf( outfile,  "imax 1  number of channels\n" );

	fprintf( outfile,  "jmax %d  number of backgrounds (", nBkgs );
	for (TString bkgName : bkgs ) fprintf( outfile, "%s, ",  bkgName.Data() );   // Eventually this will change to 4 (1l, 2l, wjets, znunu)
	fprintf( outfile,  ")\n" );

	fprintf( outfile,  "kmax %d  number of uncertainties\n", 2*nSamples-1 ); // For now, we've got statistical and systematic for each sample
	fprintf( outfile,  "---\n" );

	fprintf( outfile,  "# Now list the number of events observed (namely, zero)\n" ); // Eventually change this once I'm running over data
	fprintf( outfile,  "bin %d\n", bin );
	fprintf( outfile,  "observation 0\n" );
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



	//Now the tricky part: getting the expectation from each of the histograms!
	fprintf( outfile,  "rate     ");
	TH1D* h_yield;
	double yield;

	for( TString sample : samples ) {
	  h_yield = (TH1D*)histfile->Get( Form("sregion_%s", sample.Data()) );
	  yield   = h_yield->GetBinContent(bin);
	  fprintf( outfile,  " %10f", yield );
	}
	fprintf( outfile, "\n" );

	fprintf( outfile,  "---\n" );
	fprintf( outfile,  "# Now we list the independent sources of uncertainty (syst. error), and which samples they affect\n" ); // Can also use this bit to do stat errors
	fprintf( outfile,  "---\n" );



	// Make the table of statistical and systematic errors on each process
	for( int sampleIdx = 0; sampleIdx<nSamples; sampleIdx++ ) {

	  // Get the stat/syst errors from the appropriate histograms
	  TH1D* h_stat = (TH1D*)histfile->Get("StatUnc_"+samples.at(sampleIdx));
	  TH1D* h_syst = (TH1D*)histfile->Get("SystUnc_"+samples.at(sampleIdx));
	  double statErr = 1.0 + h_stat->GetBinContent(bin);
	  double systErr = 1.0 + h_syst->GetBinContent(bin);

	  // Generate strings for the names of these uncertainties
	  char systname[25];
	  char statname[25];
	  sprintf( systname, "%sSyst", samples.at(sampleIdx).Data() );
	  sprintf( statname, "%sStat", samples.at(sampleIdx).Data() );


	  // Print out the row for the systematic uncertainty
	  if( !(samples.at(sampleIdx).Contains("stop")
			|| samples.at(sampleIdx).Contains("signal")
			|| samples.at(sampleIdx).Contains("2tt")
			) ) {
		fprintf( outfile,   "%-18s  lnN ", systname );
		for( int j=0; j<nSamples; j++ ) {
		  if( j == sampleIdx )  fprintf( outfile, "  %4.2f  ", systErr);
		  else fprintf( outfile,  "   -    " ); // May need to adjust the padding here to match the length of the number printed above^
		}
		fprintf( outfile, "\n" );
	  }

	  // Print out the row for the statistical uncertainty
	  fprintf( outfile,   "%-18s  lnN ", statname );
	  for( int j=0; j<nSamples; j++ ) {
		if( j == sampleIdx )  fprintf( outfile, "  %4.2f  ", statErr);
		else fprintf( outfile,  "   -    " ); // May need to adjust the padding here to match the length of the number printed above^
	  }
	  fprintf( outfile, "\n" );


	} // End loop over processes (table rows)

	histfile->Close();
	delete histfile;

	fprintf( outfile,  "---\n" );

	///////////////////
	// Now restore the output to the terminal
	
	// cout.rdbuf( restore );
	// outfile.close();

	fclose(outfile);



  } // End loop over signal regions

}
