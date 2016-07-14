#include "dataMCplotMaker.h"
#include "analysis.h"
#include "sample.h"

#include "TFile.h"
// #include "TColor.h"

#include <cmath>

using namespace std;

void makeStack( analysis* myAnalysis) {

  TFile* plotfile = new TFile( myAnalysis->GetFileName(), "READ" );

  vector<TString> varNames;       vector<TString> axisLabels; //variables
  varNames.push_back("mt"     );  axisLabels.push_back("M_{T}");
  varNames.push_back("met"    );  axisLabels.push_back("MET");
  varNames.push_back("mt2w"   );  axisLabels.push_back("MT2W");
  varNames.push_back("chi2"   );  axisLabels.push_back("#chi^{2}");
  varNames.push_back("htratio");  axisLabels.push_back("H_{T} ratio");
  varNames.push_back("mindphi");  axisLabels.push_back("#Delta#phi");
  varNames.push_back("ptb1"   );  axisLabels.push_back("p_{T}");
  varNames.push_back("drlb1"  );  axisLabels.push_back("#DeltaR");
  varNames.push_back("ptlep"  );  axisLabels.push_back("p_{T}");
  varNames.push_back("metht"  );  axisLabels.push_back("MET/sqrt H_{T}");
  varNames.push_back("dphilw" );  axisLabels.push_back("#Delta#phi");
  varNames.push_back("bkgtype");  axisLabels.push_back("Background type");
  varNames.push_back("njets"  );  axisLabels.push_back("Number of jets");
  varNames.push_back("nbtags" );  axisLabels.push_back("Number of b-tags");
  varNames.push_back("ptj1"   );  axisLabels.push_back("p_{T}");
  varNames.push_back("j1btag" );  axisLabels.push_back("b-tagged?");


  // Get the signal region names
  vector<TString>  regNames = myAnalysis->GetSigRegionLabelsAll();
  vector<sigRegion> regions = myAnalysis->GetSigRegionsAll();

  // Get sample titles and colors from the "analysis" object
  vector<string> bkg_titles = myAnalysis->GetBkgNamesLegend();
  vector<string> sig_titles = myAnalysis->GetSignalNamesLegend();
  vector<short int> colors  = myAnalysis->GetColors();

  // Convert the lumi to a TString
  char tmpStr[10];
  sprintf( tmpStr, "%f", myAnalysis->GetLumi() );
  TString lumi(tmpStr);



  // Loop over all signal regions, and all the variables we're plotting
  for( uint j=0; j<regNames.size(); j++ ) {
	for( uint i=0; i<varNames.size(); i++ ) {

	  vector<TH1F*> bkgs;
	  vector<TH1F*> sigs;
	  TH1F* data;

	  // Retrieve the histograms for each background
	  for( TString sampleName : myAnalysis->GetBkgLabels() ) {
		TString plotname = varNames.at(i) + "_" + sampleName + "_" + regNames.at(j);
		TH1F*   histo    = (TH1F*)plotfile->Get(plotname);
		bkgs.push_back(histo);
	  }

	  // Retrieve the histograms for each signal
	  for( TString sampleName : myAnalysis->GetSignalLabels() ) {
		TString plotname = varNames.at(i) + "_" + sampleName + "_" + regNames.at(j);
		TH1F*   histo    = (TH1F*)plotfile->Get(plotname);
		sigs.push_back(histo);
	  }

	  // If there's a data sample, get the histogram
	  if( myAnalysis->HasData() ) {
		TString plotname = varNames.at(i) + "_" + myAnalysis->GetData()->GetLabel() + "_" + regNames.at(j);
		data = (TH1F*)plotfile->Get(plotname);
	  }
	  else data = new TH1F("", "", 1, 0, 1);

	  // Get the title and subtitle for the plot, and do a little sanity check
	  TString plotTitle;
	  TString plotSubTitle = "Region: " + regions.at(j).GetRootName();

	  if(   myAnalysis->GetNsignals() > 0 ) plotTitle = sigs.at(0)->GetTitle();
	  else if( myAnalysis->GetNbkgs() > 0 ) plotTitle = bkgs.at(0)->GetTitle();
	  else if( myAnalysis->HasData()      ) plotTitle = data->GetTitle();
	  else {
		cout << "Error in makeStack.cc: The analysis object has no samples!" << endl;
		throw(5);
	  }
	  


	  // Option string for the stack maker
	  TString optString = "--energy 13 --lumi " + lumi + " --xAxisLabel "+axisLabels.at(i)+" --xAxisUnit --outputName plots/stack_" + varNames.at(i) + "_" + regNames.at(j); // + " --png";

	  // Run the big tamale...
	  dataMCplotMaker( data,
					   bkgs,
					   bkg_titles,
					   plotTitle.Data(),
					   plotSubTitle.Data(),
					   optString.Data(),
					   sigs,
					   sig_titles,
					   colors );

	} // End loop over variables to plot
  } // End loop over signal regions


  // Also make a stack for the yields by signal region (i.e. the "final plot")
  vector<TH1F*> bkgs;
  vector<TH1F*> sigs;
  TH1F* data;
  for( TString sampleName : myAnalysis->GetBkgLabels()    ) bkgs.push_back(  (TH1F*)plotfile->Get("srYields_"+sampleName)  );
  for( TString sampleName : myAnalysis->GetSignalLabels() ) sigs.push_back(  (TH1F*)plotfile->Get("srYields_"+sampleName)  );
  if( myAnalysis->HasData() ) data = (TH1F*)plotfile->Get( "srYields_" + myAnalysis->GetData()->GetLabel() );
  else data = new TH1F("", "", 1, 0, 1);
  TString optString = "--energy 13 --lumi " + lumi + " --xAxisLabel Signal Region --xAxisUnit --outputName plots/stack_srYields";
  dataMCplotMaker( data, bkgs, bkg_titles, "Yields by signal region", "", optString.Data(), sigs, sig_titles, colors );

  plotfile->Close();
  delete plotfile;

}
