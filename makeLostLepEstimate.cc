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


  // Get lost lepton background yields from MC in signal regions
  TFile* uncFile = new TFile("uncertSR.root", "READ");
  for( uint i=0; i<srnames.size(); i++ ) {
	TH1D* histo = (TH1D*)uncFile->Get("evttype_"+srnames.at(i));
	h_srMC->SetBinContent(i+1, histo->GetBinContent(4));
	h_srMC->SetBinError(i+1, histo->GetBinError(4));
  }
  uncFile->Close();
  delete uncFile;
  outFile->cd();

  // Get total yields from MC in 2-lep control regions
  for( TString sampleName : crAnalysis->GetBkgLabels() ) {
	TH1D* histo = (TH1D*)crHistFile->Get("srYields_"+sampleName);
	h_crMC->Add( histo );
  }

  // Now do the division, M^SR / M^CR
  TH1D* h_mcRatio = (TH1D*)h_srMC->Clone("mcRatioLostLep");
  h_mcRatio->SetTitle( "SR/CR ratio by signal region" );
  for( uint i=0; i<nSRegions; i++ ) h_crMC->GetXaxis()->SetBinLabel( i+1, srnames.at(i) ); // equalize bin names
  h_mcRatio->Divide( h_crMC );


  // Multiply M/M by N^CR
  TH1D* h_crData;
  TH1D* h_bkgEstimate;

  if( crAnalysis->HasData() ) h_crData = (TH1D*)crHistFile->Get("srYields_"+crAnalysis->GetData()->GetLabel());
  else {
	h_crData = (TH1D*)h_crMC->Clone("crData");
	cout << "Warning in makeLostLepEstimate.cc: No data sample found in control region. Using MC as a dummy instead." << endl;
  }
  h_crData->SetTitle( "Control region yields from data" );

  h_bkgEstimate = (TH1D*)h_crData->Clone( "lostLepBkg" );
  h_bkgEstimate->SetTitle( "Lost lepton background estimate" );
  h_bkgEstimate->Multiply( h_mcRatio );


  // Write everything to a file
  h_mcRatio->Write();
  h_bkgEstimate->Write();
  cout << "Lost lepton background estimate saved in " << outFile->GetName() << "." << endl;

  // Clean up
  outFile->Close();
  crHistFile->Close();
  srHistFile->Close();

  delete outFile;
  delete crHistFile;
  delete srHistFile;

}
