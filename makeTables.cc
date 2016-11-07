#include <iostream>
#include <cmath>
// #include <string>

#include "TFile.h"
#include "TH1.h"

#include "analysis.h"
#include "sample.h"

using namespace std;

void makeTables( analysis* myAnalysis ) {

	TFile* infile = new TFile( myAnalysis->GetPlotFileName(), "READ" );
	if( infile->IsZombie() ) {
		cout << "Error in makeTables! Couldn't open " << infile->GetName() << "!" << endl;
		return;
	}

	vector<TString> decayNames;
	decayNames.push_back("Z $\\rightarrow \\nu\\nu$");
	decayNames.push_back("$\\geq$2 leptons");
	decayNames.push_back("1 lepton from top");
	decayNames.push_back("1 lepton from W");
	decayNames.push_back("Other");

	TString dummySampleName = myAnalysis->GetAllSamples().at(0)->GetLabel();
	TString dummyRegionName = myAnalysis->GetSigRegionLabelsAll().at(0);

	uint binOffset = 0;


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for( vector<sigRegion*> regList : myAnalysis->GetSigRegions() ) {

		double yield, error;

		// Retrieve histograms for the total bkg yields by decay type
		uint nRegions = regList.size();
		vector<TH1D*> yieldsByDecayType;
		for( uint i=0; i<nRegions; i++ ) yieldsByDecayType.push_back( (TH1D*)infile->Get("evttype_" + regList.at(i)->GetLabel()) );

		// Begin LaTeX table
		cout << "\\begin{tabular}{ l | ";
		for( sigRegion* sReg : regList ) cout << "c ";
		cout << "}" << endl;

		// Print column headers
		cout << " Sample  ";
		for( sigRegion* sReg : regList ) cout << "  &  " << sReg->GetTableName();
		cout << " \\\\ \\hline" << endl;


		/////////////////////////////////////////////////////
		// If we have a data sample, print a table row for it
		if( myAnalysis->HasData() ) {

			sample* data = myAnalysis->GetData();
			printf( "%28s  ", data->GetTableName().Data() );       // Print start of row
			TH1D* h_yields = (TH1D*)infile->Get( "srYields_" + data->GetLabel() ); // Retrieve yield histo for this sample

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
			TH1D* h_yields = (TH1D*)infile->Get( "srYields_" + thisSample->GetLabel() ); // Retrieve yield histo for this sample

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
			TH1D* h_yields = (TH1D*)infile->Get( "srYields_" + thisSample->GetLabel() ); // Retrieve yield histo for this sample

			// Read in yields and errors, and print out another cell in the table row
			for( uint i=1; i<=nRegions; i++ ) {
				yield = h_yields->GetBinContent( i + binOffset );
				error = h_yields->GetBinError( i + binOffset );
				printf( "  &   %8.3f $\\pm$ %6.3f", yield, error );
			}
			cout << "  \\\\" << endl;
		} // End loop over background samples

		binOffset += nRegions;



		/////////////////////////////////////////////////////////////////////////////////////
		// If we have background samples, then print a row for the total background yield...

		if( myAnalysis->GetNbkgs() > 0 ) {
			cout << "\\hline" << endl;

			cout << "            Total Background  ";
			for( uint i=1; i<=nRegions; i++ ) {
				yield = yieldsByDecayType.at(i-1)->IntegralAndError( 3, 6, error );
				printf( "  &   %8.3f $\\pm$ %6.3f", yield, error);
			}
			cout << "  \\\\" << endl;
			cout << "\\hline \\hline" << endl;

			/////////////////////////////////////////////////////////////////////////////////
			// ...and make the lower half of the table, with the background yields by type

			for( uint i=0; i<decayNames.size()-1; i++ ) {

				printf( "%28s  ", decayNames.at(i).Data() );

				for( uint j=0; j<nRegions; j++ ) {
					yield = yieldsByDecayType.at(j)->GetBinContent(i+3);
					error = yieldsByDecayType.at(j)->GetBinError(i+3);
					printf( "  &   %8.3f $\\pm$ %6.3f", yield, error );
				}
				cout << "  \\\\" << endl;

			} // end loop over signal regions
		} // end "if we have background samples"


		///////////////////////////////////////////////
		// End LaTeX table
		cout << "\\end{tabular}\n\n" << endl;

	} // end loop over lists of signal regions

	infile->Close();
	delete infile;

}
