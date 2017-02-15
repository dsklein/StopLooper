#include <iostream>

#include "analysis.h"
#include "sigRegion.h"
#include "systematic.h"

#include "TAxis.h"
#include "TFile.h"
#include "TH1D.h"

using namespace std;


TH1D* convertHisto( TFile* srcFile, TString oldname, TString newname, vector<TString> labels );


void makeZNuNuEstimate( analysis* znunuAnalysis ) {

	TH1::SetDefaultSumw2();

	TFile* infile = new TFile( "ZNuNu_BkgEst.root", "READ" );
	TFile* outfile = new TFile( "zNuNuEstimate.root", "RECREATE" );
	outfile->cd();

	vector<TString> regionList = znunuAnalysis->GetSigRegionLabelsAll();

	vector<TH1D*> newHistos;

	// Process all of Marketa's histograms into my style
	newHistos.push_back( convertHisto( infile, "yield", "zNuNuBkg", regionList ) );
	newHistos.push_back( convertHisto( infile, "lepSFDN", "variation_lepSFdown", regionList ) );
	newHistos.push_back( convertHisto( infile, "lepSFUP", "variation_lepSFup", regionList ) );
	newHistos.push_back( convertHisto( infile, "btagLightDN", "variation_btagLFdown", regionList ) );
	newHistos.push_back( convertHisto( infile, "btagLightUP", "variation_btagLFup", regionList ) );
	newHistos.push_back( convertHisto( infile, "btagHeavyDN", "variation_btagHFdown", regionList ) );
	newHistos.push_back( convertHisto( infile, "btagHeavyUP", "variation_btagHFup", regionList ) );
	newHistos.push_back( convertHisto( infile, "PUdown", "variation_PUup", regionList ) );
	newHistos.push_back( convertHisto( infile, "PUup", "variation_PUdown", regionList ) );
	newHistos.push_back( convertHisto( infile, "pdfDN", "variation_pdfdown", regionList ) );
	newHistos.push_back( convertHisto( infile, "pdfUP", "variation_pdfup", regionList ) );
	newHistos.push_back( convertHisto( infile, "alphaSDN", "variation_alphaSdown", regionList ) );
	newHistos.push_back( convertHisto( infile, "alphaSUP", "variation_alphaSup", regionList ) );
	newHistos.push_back( convertHisto( infile, "Q2DN", "variation_qSquareddown", regionList ) );
	newHistos.push_back( convertHisto( infile, "Q2UP", "variation_qSquaredup", regionList ) );
	newHistos.push_back( convertHisto( infile, "jesDN", "variation_JESdown", regionList ) );
	newHistos.push_back( convertHisto( infile, "jesUP", "variation_JESup", regionList ) );
	newHistos.push_back( convertHisto( infile, "normalizationDN", "variation_normdown", regionList ) );
	newHistos.push_back( convertHisto( infile, "normalizationUP", "variation_normup", regionList ) );
	newHistos.push_back( convertHisto( infile, "total", "total_uncert", regionList ) );

	// Write all the new histograms to my file
	for( TH1D* thisHist : newHistos ) {
		thisHist->Write();
		cout << "Systematic " << thisHist->GetName() << " saved in " << gFile->GetName() << "." << endl;
	}


	// Cleanup
	delete outfile;
	delete infile;
}



// Function to convert histograms from Marketa's style and naming conventions to my own
TH1D* convertHisto( TFile* srcFile, TString oldname, TString newname, vector<TString> labels ) {

	TH1D* oldhist = (TH1D*)srcFile->Get(oldname);
	if( oldhist == NULL ) {
		cout << "Error in makeZNuNuEstimate! Histogram " << oldname << " does not exist in file " << srcFile->GetName() << "!" << endl;
		throw(5);
	}

	// Check to make sure the number of bins matches
	int nBinsNew = int(labels.size());
	int nBinsOld = oldhist->GetNbinsX();
	if( nBinsOld != nBinsNew ) {
		cout << "Error in makeZNuNuEstimate! Histogram \"" << oldname << "\" has " << nBinsOld << " bins, and we have " << nBinsNew << " signal regions!" << endl;
		throw(5);
	}

	// Copy bin contents and errors to new histogram, and relabel bins
	TH1D* newhist = new TH1D( newname, newname, nBinsNew, 0.5, float(nBinsNew)+0.5 );
	TAxis* axis = newhist->GetXaxis();
	for( int i=1; i<=nBinsNew; i++ ) {
		newhist->SetBinContent( i, oldhist->GetBinContent(i) );
		newhist->SetBinError( i, oldhist->GetBinError(i) );
		axis->SetBinLabel( i, labels.at(i-1) );
	}

	return newhist;
}
