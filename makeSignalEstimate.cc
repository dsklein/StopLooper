#include <iostream>

#include "analysis.h"
#include "sample.h"
#include "sigRegion.h"
#include "systematic.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"

using namespace std;


void makeSignalEstimate( analysis* sigAnalysis ) {

	TH1::SetDefaultSumw2();

	// Open input files and output file
	TFile* sigHistFile = new TFile( sigAnalysis->GetPlotFileName(), "READ" );
	TFile* sigSystFile = new TFile( sigAnalysis->GetSystFileName(), "READ" );
	TFile* sigJesFile  = new TFile( "jes_sr.root", "READ" );
	TFile* lostlepFile = new TFile( "lostlepEstimates.root", "READ" );
	TFile* onelepwFile = new TFile( "onelepwEstimates.root", "READ" );

	// Check for bad files
	for( TFile* thisFile : {sigHistFile,sigSystFile,sigJesFile,lostlepFile,onelepwFile} ) {
		if( thisFile->IsZombie() ) {
			cout << "Error in makeSignalEstimate! Couldn't open file " << thisFile->GetName() << "!" << endl;
			return;
		}
	}

	TFile* outFile = new TFile( "signalEstimates.root", "RECREATE" );
	outFile->cd();

	// Define some useful variables
	vector<TString> srnames = sigAnalysis->GetSigRegionLabelsAll();
	vector<vector<sigRegion*> > sigRegionList = sigAnalysis->GetSigRegions();


	// Loop over the signal regions
	for( TString regName : sigAnalysis->GetSigRegionLabelsAll() ) {

		// Copy the nominal yields to the output file
		TH2D* h_sigYield = (TH2D*)sigHistFile->Get("sigyields_"+regName)->Clone("nominal_"+regName);
		h_sigYield->Write();

		// Copy the varied yields to the output file
		for( systematic* thisSys : sigAnalysis->GetSystematics(true) ) {
			TString suffix = "_" + thisSys->GetNameLong();
			TFile* infile = suffix.Contains("JES") ? sigJesFile : sigSystFile;
			TH2D* h_tmp = (TH2D*)infile->Get("sigyields_"+regName+suffix)->Clone("variation_"+regName+suffix);
			h_tmp->Write();
		}

		// Subtract signal contamination from lostlep and 1lW backgrounds
		TH2D* h_contam_subtracted = (TH2D*)h_sigYield->Clone("contamSubtracted_"+regName);
		TH2D* h_lostlep_contam    = (TH2D*)lostlepFile->Get("sigContam_"+regName);
		TH2D* h_onelepw_contam    = (TH2D*)onelepwFile->Get("sigContam_"+regName);
		h_contam_subtracted->Add( h_lostlep_contam, -1. );
		h_contam_subtracted->Add( h_onelepw_contam, -1. );
		h_contam_subtracted->Write();

		// Also do the fastsim genMET averaging thing


		cout << "Signal yields and variations for region " << regName << " stored in file " << gFile->GetName() << "." << endl;
	} // End loop over signal regions

	// Clean up
	outFile->Close();
	sigSystFile->Close();
	sigHistFile->Close();
	sigJesFile->Close();

	delete outFile;
	delete sigSystFile;
	delete sigHistFile;
	delete sigJesFile;

}
