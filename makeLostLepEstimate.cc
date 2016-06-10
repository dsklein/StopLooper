#include <iostream>

#include "analysis.h"
#include "sample.h"

#include "TFile.h"
#include "TH1.h"

using namespace std;

void makeLostLepEstimate( analysis* srAnalysis, analysis* crAnalysis ) {

  TH1::SetDefaultSumw2();

  vector<TString> srnames = srAnalysis->GetSigRegionsAll();
  vector<TString> crnames = crAnalysis->GetSigRegionsAll();
  uint nSRegions = srnames.size();
  uint nCRegions = crnames.size();
  if( nSRegions != nCRegions ) {
	cout << "Error in makeLostLepEstimate: Different number of signal and control regions!" << endl;
	return;
  }


  // Open input files and output file
  TFile* srHistFile = new TFile( srAnalysis->GetFileName(), "READ" );
  TFile* crHistFile = new TFile( crAnalysis->GetFileName(), "READ" );
  TFile* outFile    = new TFile( "bkgEstimates.root",   "RECREATE" );


  // Declare histograms that will hold the SR/CR yields
  TH1D* h_srMC   = new TH1D( "srMC"  , "Signal region yields from MC"   , nSRegions, 0.5, float(nSRegions)+0.5 );
  TH1D* h_crMC   = new TH1D( "crMC"  , "Control region yields from MC"  , nCRegions, 0.5, float(nCRegions)+0.5 );
  for( uint i=0; i<nSRegions; i++ ) h_srMC->GetXaxis()->SetBinLabel( i+1, srnames.at(i) );
  for( uint i=0; i<nCRegions; i++ ) h_crMC->GetXaxis()->SetBinLabel( i+1, crnames.at(i) );


  // Get total background yields from MC in signal regions
  for( TString sampleName : srAnalysis->GetBkgLabels() ) {
	TH1D* histo = (TH1D*)srHistFile->Get("srYields_"+sampleName);
	h_srMC->Add( histo );
  }

  // Get total background yields from MC in control regions
  for( TString sampleName : crAnalysis->GetBkgLabels() ) {
	TH1D* histo = (TH1D*)crHistFile->Get("srYields_"+sampleName);
	h_crMC->Add( histo );
  }

  // Now do the division, M^SR / M^CR
  TH1D* h_mcRatio = (TH1D*)h_srMC->Clone("mcRatioLostLep");
  h_mcRatio->SetTitle( "SR/CR ratio by signal region" );
  for( uint i=0; i<nSRegions; i++ ) h_crMC->GetXaxis()->SetBinLabel( i+1, srnames.at(i) ); // equalize bin names
  h_mcRatio->Divide( h_crMC );


  // Set up placeholder code to multiply M/M by N^CR
  //   Once you know how, properly account for uncertainties

  TH1D* h_crData;
  TH1D* h_bkgEstimate;

  if( crAnalysis->HasData() ) h_crData = (TH1D*)crHistFile->Get("srYields_"+crAnalysis->GetData()->GetLabel());
  else {
	h_crData = (TH1D*)h_crMC->Clone("crData");
	cout << "Warning in makeLostLepEstimate.cc: No data sample found in control region. Using MC as a dummy instead." << endl;
  }
  h_crData->SetTitle( "Control region yields from data" );

  h_bkgEstimate = (TH1D*)h_crMC->Clone( "lostLepBkg" );
  h_bkgEstimate->SetTitle( "Lost lepton background estimate" );
  h_bkgEstimate->Multiply( h_mcRatio );


  // Write everything to a file
  h_mcRatio->Write();
  h_bkgEstimate->Write();

  // Clean up
  outFile->Close();
  crHistFile->Close();
  srHistFile->Close();

  delete outFile;
  delete crHistFile;
  delete srHistFile;

}
