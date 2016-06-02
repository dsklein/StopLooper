#include <iostream>
#include <cmath>
// #include <string>

#include "TFile.h"
#include "TH1.h"

#include "analysis.h"
#include "sample.h"

using namespace std;

void makeTables( analysis* myAnalysis ) {

  TFile* infile = new TFile( myAnalysis->GetFileName(), "READ" );

  vector<TString> decayNames;
  decayNames.push_back("1 lepton");
  decayNames.push_back("$\\geq$2 leptons");
  decayNames.push_back("Z $\\rightarrow \\nu\\nu$");
  decayNames.push_back("Other");

  TString dummySampleName = myAnalysis->GetBkgLabels().at(0);
  TString dummyRegionName = myAnalysis->GetSigRegionsAll().at(0);

  // Keep a histogram for the total yields
  TH1D* h_totals_sregion = (TH1D*)infile->Get("sregion_" + dummySampleName)->Clone("signal_region_totals");

  uint binOffset = 0;


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for( vector<TString> regNameList : myAnalysis->GetSigRegions() ) {

	double yield, error;
	h_totals_sregion->Reset();

	// Prep histograms that will store the total bkg yields by decay type
	uint nRegions = regNameList.size();
	TH1D* yieldsByDecayType[nRegions];
	for( uint i=0; i<nRegions; i++ ) {
	  yieldsByDecayType[i] = (TH1D*)infile->Get("bkgtype_" + dummySampleName + "_" + dummyRegionName)->Clone( "decayType_" + regNameList.at(i) );
	  yieldsByDecayType[i]->Reset();
	}

	// Begin LaTeX table
	cout << "\\begin{tabular}{ l | ";
	for( TString srName : regNameList ) cout << "c ";
	cout << "}" << endl;

	// Print column headers
	cout << " Sample  ";
	for( TString srName : regNameList ) cout << "  &  " << srName;
	cout << " \\\\ \\hline" << endl;


	/////////////////////////////////////////////////////
	// If we have a data sample, print a table row for it
	if( myAnalysis->HasData() ) {

	  sample* data = myAnalysis->GetData();
	  printf( "%28s  ", data->GetTableName().Data() );       // Print start of row
	  TH1D* h_yields = (TH1D*)infile->Get( "sregion_" + data->GetLabel() ); // Retrieve yield histo for this sample

	  // Read in yields and errors, and print out another cell in the table row
	  for( uint i=1; i<=nRegions; i++ ) {
		yield = h_yields->GetBinContent( i + binOffset );
		error = h_yields->GetBinError( i + binOffset );
		printf( "  &   %8.3f $\\pm$ %6.3f", yield, error );
	  }
	  cout << "  \\\\" << endl;
	  cout << "\\hline" << endl;
	}

	///////////////////////////////////////////////////////////
	// Loop over signal samples, and print a table row for each
	for( sample* thisSample : myAnalysis->GetSignals() ) {

	  printf( "%28s  ", thisSample->GetTableName().Data() );       // Print start of row
	  TH1D* h_yields = (TH1D*)infile->Get( "sregion_" + thisSample->GetLabel() ); // Retrieve yield histo for this sample

	  // Read in yields and errors, and print out another cell in the table row
	  for( uint i=1; i<=nRegions; i++ ) {
		yield = h_yields->GetBinContent( i + binOffset );
		error = h_yields->GetBinError( i + binOffset );
		printf( "  &   %8.3f $\\pm$ %6.3f", yield, error );
	  }
	  cout << "  \\\\" << endl;
	} // End loop over signal samples


	if( myAnalysis->GetNsignals() > 0 ) cout << "\\hline" << endl;

	///////////////////////////////////////////////////////////////
	// Loop over background samples, and print a table row for each
	for( sample* thisSample : myAnalysis->GetBkgs() ) {

	  printf( "%28s  ", thisSample->GetTableName().Data() );       // Print start of row
	  TH1D* h_yields = (TH1D*)infile->Get( "sregion_" + thisSample->GetLabel() ); // Retrieve yield histo for this sample
	  h_totals_sregion->Add( h_yields );

	  // Read in yields and errors, and print out another cell in the table row
	  // Also add yields by background type to running tally
	  for( uint i=1; i<=nRegions; i++ ) {
		yield = h_yields->GetBinContent( i + binOffset );
		error = h_yields->GetBinError( i + binOffset );
		printf( "  &   %8.3f $\\pm$ %6.3f", yield, error );
		
		TH1D* h_decayType = (TH1D*)infile->Get( "bkgtype_" + thisSample->GetLabel() + "_" + regNameList.at(i-1) );
		yieldsByDecayType[i-1]->Add(h_decayType);
	  }
	  cout << "  \\\\" << endl;
	} // End loop over background samples



	cout << "\\hline" << endl;

	////////////////////////////////////////
	// Print row for total background yield

	cout << "            Total Background  ";
	for( uint i=1; i<=nRegions; i++ ) {
	  yield = h_totals_sregion->GetBinContent(i + binOffset);
	  error = h_totals_sregion->GetBinError(i + binOffset);
	  printf( "  &   %8.3f $\\pm$ %6.3f", yield, error);
	}
	cout << "  \\\\" << endl;
	cout << "\\hline \\hline" << endl;

	binOffset += nRegions;

	/////////////////////////////////////////////////////////////////////////////////
	// Make the lower half of the table, with the rows = final states

	for( uint i=0; i<decayNames.size()-1; i++ ) {

	  printf( "%28s  ", decayNames.at(i).Data() );

	  for( uint j=0; j<nRegions; j++ ) {
		yield = yieldsByDecayType[j]->GetBinContent(i+1);
		error = yieldsByDecayType[j]->GetBinError(i+1);
		printf( "  &   %8.3f $\\pm$ %6.3f", yield, error );
	  }
	  cout << "  \\\\" << endl;

	} // end loop over signal regions

	cout << "\\end{tabular}\n\n" << endl;


  } // end loop over lists of signal regions
 
  infile->Close();
  delete infile;

}
