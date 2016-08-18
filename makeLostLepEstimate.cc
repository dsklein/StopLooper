#include <iostream>

#include "analysis.h"
#include "sample.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"

using namespace std;

void makeLostLepEstimate( analysis* srAnalysis, analysis* crAnalysis ) {

  TH1::SetDefaultSumw2();

  vector<TString> srnames = srAnalysis->GetSigRegionLabelsAll();
  vector<TString> crnames = crAnalysis->GetSigRegionLabelsAll();
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
  for( uint i=0; i<srnames.size(); i++ ) {
	TH1D* histo = (TH1D*)srHistFile->Get("evttype_"+srnames.at(i));
	h_srMC->SetBinContent(i+1, histo->GetBinContent(4));
	h_srMC->SetBinError(i+1, histo->GetBinError(4));
  }

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


  // Get data yields in CRs, and multiply by M/M
  TH1D* h_crData;
  TH1D* h_bkgEstimate;

  if( crAnalysis->HasData() ) h_crData = (TH1D*)crHistFile->Get("srYields_"+crAnalysis->GetData()->GetLabel())->Clone("crData");
  else {
	h_crData = (TH1D*)h_crMC->Clone("crData");
	cout << "\nWarning in makeLostLepEstimate.cc: No data sample found in control region. Using MC as a dummy instead." << endl;
  }
  h_crData->SetTitle( "Control region yields from data" );
  h_crData->Write();
  for( uint i=0; i<nSRegions; i++ ) h_crData->GetXaxis()->SetBinLabel( i+1, srnames.at(i) ); // equalize bin names

  h_bkgEstimate = (TH1D*)h_crData->Clone( "lostLepBkg" );
  h_bkgEstimate->SetTitle( "Lost lepton background estimate" );
  h_bkgEstimate->Multiply( h_mcRatio );

  // Write everything to a file
  h_mcRatio->Write();
  h_bkgEstimate->Write();
  cout << "Lost lepton background estimate saved in " << outFile->GetName() << "." << endl;


  ////////////////////////////////////////////////////////////////
  // Now do some calculations for the systematics...

  // Calculate the signal contamination in the CRs
  if( crAnalysis->GetNsignals() >= 1 ) {
	for( uint i=0; i<nCRegions; i++ ) {
	  double ratio = h_mcRatio->GetBinContent( i+1 );
	  double ratio_err = h_mcRatio->GetBinError( i+1 );
	  TH2D* h_contam = (TH2D*)crHistFile->Get("sigyields_"+crnames.at(i))->Clone("sigContam_"+srnames.at(i));
	  for( int bin=0; bin<h_contam->GetNcells(); bin++ ) {
		double yield = h_contam->GetBinContent(bin);
		double error = h_contam->GetBinError(bin);
		if( yield < 0.000001 ) continue;
		h_contam->SetBinContent( bin, yield*ratio );
		h_contam->SetBinError( bin, sqrt( yield*yield*ratio_err*ratio_err + ratio*ratio*error*error ) ); // Gaussian error propagation
	  }
	  h_contam->Write();
	}
	cout << "Signal contamination estimate saved in " << outFile->GetName() << "." << endl;
  }

  // Isolate the uncertainties due to signal stats and MC stats
  TH1D* h_datastats = (TH1D*)h_mcRatio->Clone("estimate_datastats");
  TH1D* h_mcstats   = (TH1D*)h_crData->Clone("estimate_mcstats");
  h_datastats->SetTitle( "Background estimate with uncertainty from data stats only" );
  h_mcstats->SetTitle( "Background estimate with uncertainty from MC stats only" );
  for( uint i=1; i<=nSRegions; i++ ) {
  	h_datastats->SetBinError(i, 0.);
  	h_mcstats->SetBinError(i, 0.);
  }
  h_datastats->Multiply( h_crData );
  h_mcstats->Multiply( h_mcRatio );
  h_datastats->Write();
  h_mcstats->Write();


  ////////////////////////////////////////////////////////////////////
  // Print out a table with the details of the lost lepton estimate

  //  Print table header
  cout << "\nLost lepton background estimate (stat errors only)\n" << endl;
  cout << "\\begin{tabular}{ | l | c | c | c | }" << endl;
  cout << "\\hline" << endl;
  cout << "Signal region  &  $N^{CR}$  &  Transfer factor  &  $N_{est}^{SR}$ \\\\" << endl;
  cout << "\\hline" << endl;

  // Loop through signal regions and print out table rows
  uint binOffset = 1;
  for( vector<sigRegion> sigRegs : srAnalysis->GetSigRegions() ) {
	for( uint i=0; i<sigRegs.size(); i++ ) {

	  printf( "%30s & %3d $\\pm$ %5.3f &  %5.3f $\\pm$ %5.3f  &  %5.2f $\\pm$ %5.2f  \\\\\n", sigRegs.at(i).GetTableName().Data(),
			  int(h_crData->GetBinContent(i+binOffset)), h_crData->GetBinError(i+binOffset), h_mcRatio->GetBinContent(i+binOffset),
			  h_mcRatio->GetBinError(i+binOffset), h_bkgEstimate->GetBinContent(i+binOffset), h_bkgEstimate->GetBinError(i+binOffset) );

	}
	binOffset += sigRegs.size();
	cout << "\\hline" << endl;
  }
  cout << "\\end{tabular}\n" << endl;


  // Clean up
  outFile->Close();
  crHistFile->Close();
  srHistFile->Close();

  delete outFile;
  delete crHistFile;
  delete srHistFile;

}
