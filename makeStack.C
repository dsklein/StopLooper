#include "dataMCplotMaker.h"
#include "TFile.h"
// #include "TColor.h"

#include <cmath>

using namespace std;

void makeStack() {

  TFile* plotfile = new TFile("plots.root", "READ");

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
  
  vector<TString> regNames; //signal regions
  regNames.push_back("low250");
  regNames.push_back("low300");
  regNames.push_back("low350");
  regNames.push_back("low400");
  regNames.push_back("high250");
  regNames.push_back("high300");
  regNames.push_back("high350");
  regNames.push_back("high400");
  regNames.push_back("high500");

  // Loop over all the variables we're plotting
  for( uint j=0; j<regNames.size(); j++ ) {
	for( uint i=0; i<varNames.size(); i++ ) {

	  TString name_tt2l		= varNames.at(i) + "_tt2l_"     + regNames.at(j);
	  TString name_tt1l		= varNames.at(i) + "_tt1l_"     + regNames.at(j);
	  // TString name_wb		= varNames.at(i) + "_Wb_"       + regNames.at(j);
	  // TString name_wlight	= varNames.at(i) + "_Wucsd_"    + regNames.at(j);
	  TString name_wjets    = varNames.at(i) + "_wjets_"    + regNames.at(j);
	  TString name_dy		= varNames.at(i) + "_dy_"       + regNames.at(j);
	  // TString name_stst		= varNames.at(i) + "_STstchan_" + regNames.at(j);
	  // TString name_sttw		= varNames.at(i) + "_STtWchan_" + regNames.at(j);
	  TString name_singletop= varNames.at(i) + "_singletop_" + regNames.at(j);
	  // TString name_ttw		= varNames.at(i) + "_ttw_"      + regNames.at(j);
	  // TString name_ttz		= varNames.at(i) + "_ttz_"      + regNames.at(j);
	  // TString name_vv		= varNames.at(i) + "_vv_"       + regNames.at(j);
	  TString name_rare		= varNames.at(i) + "_rare_"     + regNames.at(j);

	  TString name_stop850	= varNames.at(i) + "_stop850_"  + regNames.at(j);
	  TString name_stop650	= varNames.at(i) + "_stop650_"  + regNames.at(j);
	  TString name_stop500	= varNames.at(i) + "_stop500_"  + regNames.at(j);
	  TString name_stop425	= varNames.at(i) + "_stop425_"  + regNames.at(j);

	  // Extract a histogram for each of the backgrounds and signals
	  TH1F* h_tt2l		= (TH1F*)plotfile->Get(name_tt2l);
	  TH1F* h_tt1l		= (TH1F*)plotfile->Get(name_tt1l);
	  // TH1F* h_wb		= (TH1F*)plotfile->Get(name_wb);
	  // TH1F* h_wlight	= (TH1F*)plotfile->Get(name_wlight);
	  TH1F* h_wjets 	= (TH1F*)plotfile->Get(name_wjets);
	  TH1F* h_dy		= (TH1F*)plotfile->Get(name_dy);
	  // TH1F* h_stst		= (TH1F*)plotfile->Get(name_stst);
	  // TH1F* h_sttw		= (TH1F*)plotfile->Get(name_sttw);
	  TH1F* h_singletop	= (TH1F*)plotfile->Get(name_singletop);
	  // TH1F* h_ttw		= (TH1F*)plotfile->Get(name_ttw);
	  // TH1F* h_ttz		= (TH1F*)plotfile->Get(name_ttz);
	  // TH1F* h_vv		= (TH1F*)plotfile->Get(name_vv);
	  TH1F* h_rare		= (TH1F*)plotfile->Get(name_rare);


	  TH1F* h_stop850	= (TH1F*)plotfile->Get(name_stop850);
	  TH1F* h_stop650	= (TH1F*)plotfile->Get(name_stop650);
	  TH1F* h_stop500	= (TH1F*)plotfile->Get(name_stop500);
	  TH1F* h_stop425	= (TH1F*)plotfile->Get(name_stop425);

	  TH1F* nullData = new TH1F("", "", 1, 0, 1);

	  // String histos into vectors for Alex's plotmaking code
	  vector<TH1F*> bkgs;
	  bkgs.push_back(h_tt2l);
	  bkgs.push_back(h_tt1l);
	  // bkgs.push_back(h_wb);
	  // bkgs.push_back(h_wlight);
	  bkgs.push_back(h_wjets);
	  bkgs.push_back(h_dy);
	  // bkgs.push_back(h_stst);
	  // bkgs.push_back(h_sttw);
	  bkgs.push_back(h_singletop);
	  // bkgs.push_back(h_ttw);
	  // bkgs.push_back(h_ttz);
	  // bkgs.push_back(h_vv);
	  bkgs.push_back(h_rare);

	  vector<TH1F*> sigs;
	  sigs.push_back(h_stop850);
	  sigs.push_back(h_stop650);
	  sigs.push_back(h_stop500);
	  sigs.push_back(h_stop425);

	  // Make vectors of titles for backgrounds and signals
	  vector<string> bkg_titles;
	  bkg_titles.push_back("ttbar 2l");
	  bkg_titles.push_back("ttbar 1l");
	  // bkg_titles.push_back("W+b");
	  // bkg_titles.push_back("W+light");
	  bkg_titles.push_back("W+jets");
	  bkg_titles.push_back("Drell-Yan");
	  // bkg_titles.push_back("Single top, s/t-channels");
	  // bkg_titles.push_back("Single top, tW channel");
	  bkg_titles.push_back("Single top");
	  // bkg_titles.push_back("ttW");
	  // bkg_titles.push_back("ttZ");
	  // bkg_titles.push_back("Diboson");
	  bkg_titles.push_back("Rare");

	  vector<string> sig_titles;
	  sig_titles.push_back("stop-850-100");
	  sig_titles.push_back("stop-650-325");
	  sig_titles.push_back("stop-500-325");
	  sig_titles.push_back("stop-425-325");

	  // Make a vector of colors for the backgrounds and signals
	  vector<short int> colors;
	  colors.push_back(kGreen-4);
	  colors.push_back(kCyan-3);
	  colors.push_back(kRed-7);
	  colors.push_back(kOrange-2);
	  colors.push_back(kMagenta-5);
	  colors.push_back(kRed+2);
	  colors.push_back(kBlue+3);
	  colors.push_back(kGreen+3);
	  colors.push_back(kMagenta+3);
	  colors.push_back(kOrange+7);
	  // Need some more colors!

	  // Get the title and subtitle for the plot
	  TString plotTitle = bkgs.at(0)->GetTitle();
	  TString plotSubTitle = "Region: " + regNames.at(j);


	  // Make the options string for each stack
	  TString optString = "--energy 13 --lumi 5 --xAxisLabel "+axisLabels.at(i)+" --xAxisUnit --outputName plots/stack_" + varNames.at(i) + "_" + regNames.at(j); // + " --png";
	  // (try also --isLinear);  --legendTextSize 0.022

	  // Run the big tamale...
	  dataMCplotMaker( nullData,
					   bkgs,
					   bkg_titles,
					   plotTitle.Data(), //title
					   plotSubTitle.Data(), //subtitle
					   optString.Data(), //options
					   sigs,
					   sig_titles,
					   colors );

	} // End loop over variables to plot
  } // End loop over signal regions




}
