#include "dataMCplotMaker.h"
#include "TFile.h"
// #include "TColor.h"

#include <cmath>

using namespace std;

void makeStack() {

  TFile* plotfile = new TFile("plots.root", "READ");

  vector<TString> varNames; //variables
  varNames.push_back("mt"     );
  varNames.push_back("met"    );
  varNames.push_back("mt2w"   );
  varNames.push_back("chi2"   );
  varNames.push_back("htratio");
  varNames.push_back("mindphi");
  varNames.push_back("ptb1"   );
  varNames.push_back("drlb1"  );
  varNames.push_back("ptlep"  );
  varNames.push_back("metht"  );
  varNames.push_back("dphilw" );

  vector<TString> regNames; //signal regions
  regNames.push_back("0_150");
  regNames.push_back("0_200");
  regNames.push_back("0_250");
  regNames.push_back("0_300");
  regNames.push_back("200_150");
  regNames.push_back("200_200");
  regNames.push_back("200_250");
  regNames.push_back("200_300");

  // Loop over all the variables we're plotting
  for( uint j=0; j<regNames.size(); j++ ) {
	for( uint i=0; i<varNames.size(); i++ ) {

	  TString name_t1l   = varNames.at(i) + "_top1l_" + regNames.at(j);
	  TString name_tt2l  = varNames.at(i) + "_tt2l_"  + regNames.at(j);
	  TString name_tt1l  = varNames.at(i) + "_tt1l_"  + regNames.at(j);
	  TString name_wjets = varNames.at(i) + "_wjets_" + regNames.at(j);
	  TString name_rare  = varNames.at(i) + "_rare_"  + regNames.at(j);

	  TString name_stop850 = varNames.at(i) + "_stop850_" + regNames.at(j);
	  TString name_stop650 = varNames.at(i) + "_stop650_" + regNames.at(j);
	  TString name_stop500 = varNames.at(i) + "_stop500_" + regNames.at(j);
	  TString name_stop425 = varNames.at(i) + "_stop425_" + regNames.at(j);

	  // Extract a histogram for each of the backgrounds and signals
	  TH1F* h_tsl   = (TH1F*)plotfile->Get(name_t1l);
	  TH1F* h_ttdl  = (TH1F*)plotfile->Get(name_tt2l);
	  TH1F* h_ttsl  = (TH1F*)plotfile->Get(name_tt1l);
	  TH1F* h_wjets = (TH1F*)plotfile->Get(name_wjets);
	  TH1F* h_rare  = (TH1F*)plotfile->Get(name_rare);

	  TH1F* h_stop850 = (TH1F*)plotfile->Get(name_stop850);
	  TH1F* h_stop650 = (TH1F*)plotfile->Get(name_stop650);
	  TH1F* h_stop500 = (TH1F*)plotfile->Get(name_stop500);
	  TH1F* h_stop425 = (TH1F*)plotfile->Get(name_stop425);

	  TH1F* nullData = new TH1F("", "", 1, 0, 1);

	  // String histos into vectors for Alex's plotmaking code
	  vector<TH1F*> bkgs;
	  bkgs.push_back(h_tsl);
	  bkgs.push_back(h_ttdl);
	  bkgs.push_back(h_ttsl);
	  bkgs.push_back(h_wjets);
	  bkgs.push_back(h_rare);

	  vector<TH1F*> sigs;
	  sigs.push_back(h_stop850);
	  sigs.push_back(h_stop650);
	  sigs.push_back(h_stop500);
	  sigs.push_back(h_stop425);

	  // Make vectors of titles for backgrounds and signals
	  vector<string> bkg_titles;
	  bkg_titles.push_back("single top");
	  bkg_titles.push_back("ttbar 2l");
	  bkg_titles.push_back("ttbar 1l");
	  bkg_titles.push_back("wjets");
	  bkg_titles.push_back("rare");

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

	  // Get the title and subtitle for the plot
	  TString plotTitle = h_stop850->GetTitle();
	  TString plotSubTitle = "MT2W, MET > " + (regNames.at(j).Copy()).ReplaceAll("_", ", ");


	  // Make the options string for each stack
	  TString optString = "--energy 13 --lumi 10 --outputName plots/stack_" + varNames.at(i) + "_" + regNames.at(j); // + " --png";
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

	}
  }

  TH1D* h_yield;

  double totals[9] = {0};
  double errSqr[9] = {0};

  cout << "\n" << endl;

  // MT2W > 0 (low Delta M) region

  // Background samples
  cout << "Sample   &  MET $>$ 150  &  MET $>$ 200  &   MET $>$ 250  & MET $>$ 300  \\\\" << endl;
  cout << "\\hline" << endl;
  cout << "MT2W $>$ 0 & & & & \\\\" << endl;

  for( TString sName : {"top1l", "tt2l", "tt1l", "wjets", "rare"} ) {
	  h_yield = (TH1D*)plotfile->Get("yield_" + sName);
	  cout << sName.Data();

	  for( int j=2; j<6; j++ ) {
		double yield = h_yield->GetBinContent(j);
		double yError = h_yield->GetBinError(j);
		printf( "  &  %8.3f $\\pm$ %6.3f  ", yield, yError );
		totals[j-1] += yield;
		errSqr[j-1] += yError*yError;
	  }
	  cout << "\\\\" << endl;
  }

  // Total background
  cout << "\\hline\nTotal bkg  ";
  for( int j=1; j<5; j++ ) {
	printf( "&  %8.3f $\\pm$ %6.3f  ", totals[j], sqrt(errSqr[j]) );
  }
  cout << "\\\\\n\\hline" << endl;

  // Signal samples
  for( TString sName : {"stop850", "stop650", "stop500", "stop425"} ) {
	  h_yield = (TH1D*)plotfile->Get("yield_" + sName);
	  cout << sName.Data();

	  for( int j=2; j<6; j++ ) {
		double yield = h_yield->GetBinContent(j);
		double yError = h_yield->GetBinError(j);
		printf( "  &  %8.3f $\\pm$ %6.3f  ", yield, yError );
		totals[j-1] += yield;
		errSqr[j-1] += yError*yError;
	  }
	  cout << "\\\\" << endl;
  }










  cout << "\n" << endl;
  // MT2W > 200 (high Delta M) region

  // Background samples
  cout << "Sample   &  MET $>$ 150  &  MET $>$ 200  &   MET $>$ 250  & MET $>$ 300  \\\\" << endl;
  cout << "\\hline" << endl;
  cout << "MT2W $>$ 200 & & & & \\\\" << endl;

  for( TString sName : {"top1l", "tt2l", "tt1l", "wjets", "rare"} ) {
	  h_yield = (TH1D*)plotfile->Get("yield_" + sName);
	  cout << sName.Data();

	  for( int j=6; j<10; j++ ) {
		double yield = h_yield->GetBinContent(j);
		double yError = h_yield->GetBinError(j);
		printf( "  &  %8.3f $\\pm$ %6.3f  ", yield, yError );
		totals[j-1] += yield;
		errSqr[j-1] += yError*yError;
	  }
	  cout << "\\\\" << endl;
  }

  // Total background
  cout << "\\hline\nTotal bkg  ";
  for( int j=5; j<9; j++ ) {
	printf( "&  %8.3f $\\pm$ %6.3f  ", totals[j], sqrt(errSqr[j]) );
  }
  cout << "\\\\\n\\hline" << endl;

  // Signal samples
  for( TString sName : {"stop850", "stop650", "stop500", "stop425"} ) {
	  h_yield = (TH1D*)plotfile->Get("yield_" + sName);
	  cout << sName.Data();

	  for( int j=6; j<10; j++ ) {
		double yield = h_yield->GetBinContent(j);
		double yError = h_yield->GetBinError(j);
		printf( "  &  %8.3f $\\pm$ %6.3f  ", yield, yError );
		totals[j-1] += yield;
		errSqr[j-1] += yError*yError;
	  }
	  cout << "\\\\" << endl;
  }

}
