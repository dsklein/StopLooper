#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH1.h"

void makeTables() {

  TFile* infile = new TFile("plots.root", "READ");

  vector<TString> sampleNames;
  sampleNames.push_back("stop850");
  sampleNames.push_back("stop650");
  sampleNames.push_back("stop500");
  sampleNames.push_back("stop425");
  sampleNames.push_back("tt2l");
  sampleNames.push_back("tt1l");
  sampleNames.push_back("Wb");
  sampleNames.push_back("Wucsd");
  sampleNames.push_back("dy");
  sampleNames.push_back("ttw");
  sampleNames.push_back("ttz");
  sampleNames.push_back("STstchan");
  sampleNames.push_back("STtWchan");
  sampleNames.push_back("vv");

  vector<TString> printNames;
  printNames.push_back("T2tt (850, 100)");
  printNames.push_back("T2tt (650, 325)");
  printNames.push_back("T2tt (500, 325)");
  printNames.push_back("T2tt (425, 325)");
  printNames.push_back("$t\\bar{t} \\rightarrow 2\\ell$");
  printNames.push_back("$t\\bar{t} \\rightarrow 1\\ell$");
  printNames.push_back("W+b");
  printNames.push_back("W+light");
  printNames.push_back("Drell-Yan");
  printNames.push_back("ttW");
  printNames.push_back("ttZ");
  printNames.push_back("Single top, s/t-channels");
  printNames.push_back("Single top, tW channel");
  printNames.push_back("Diboson");

  vector<TString> regionNames;
  regionNames.push_back("low250");
  regionNames.push_back("low300");
  regionNames.push_back("low350");
  regionNames.push_back("low400");
  regionNames.push_back("high250");
  regionNames.push_back("high300");
  regionNames.push_back("high350");
  regionNames.push_back("high400");
  regionNames.push_back("high500");

  vector<TString> bkgNames;
  bkgNames.push_back("1 lepton");
  bkgNames.push_back("$\\geq$2 leptons");
  bkgNames.push_back("Z $\\rightarrow \\nu\\nu$");
  bkgNames.push_back("Other");

  // Keep a histogram for the total yields
  TH1D* h_totals_sregion = (TH1D*)infile->Get("sregion_" + sampleNames.at(0))->Clone("signal_region_totals");
  h_totals_sregion->Reset();

  // Also load in the histograms that store the decay channels
  // Loop 'i' is over signal regions; loop 'j' is over backgrounds
  TH1D* h_bgtype[9];
  for( int i=0; i<9; i++ ) {
	h_bgtype[i] = (TH1D*)infile->Get( "bkgtype_" + sampleNames.at(0) + "_" + regionNames.at(i) )->Clone("bgtype_"+regionNames.at(i));
	h_bgtype[i]->Reset();
	for( uint j=4; j<sampleNames.size(); j++ ) h_bgtype[i]->Add( (TH1D*)infile->Get("bkgtype_"+sampleNames.at(j)+"_"+regionNames.at(i)) );
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Make the first table, for the low Delta-M region

  // Make the upper part of the table, with rows = different samples, columns = signal regions

  double yield_250, err_250;
  double yield_300, err_300;
  double yield_350, err_350;
  double yield_400, err_400;

  cout << "\\begin{tabular}{ l | c c c c }" << endl;
  cout << "Low $\\Delta M$   &  250 $\\leq E_T^{miss} <$ 300  &   300 $\\leq E_T^{miss} <$ 350  &   350 $\\leq E_T^{miss} <$ 400  & $E_T^{miss} \\geq$ 400  \\\\ \\hline" << endl;

  for( uint i=0; i<sampleNames.size(); i++ ) {
	// Read in yields and errors from stored histograms
	TH1D* h_sigRegion = (TH1D*)infile->Get("sregion_" + sampleNames.at(i));
	if( !sampleNames.at(i).Contains("stop") ) h_totals_sregion->Add(h_sigRegion); // If it's not a signal sample, add it to the background totals

	yield_250 = h_sigRegion->GetBinContent( 1 );
	yield_300 = h_sigRegion->GetBinContent( 2 );
	yield_350 = h_sigRegion->GetBinContent( 3 );
	yield_400 = h_sigRegion->GetBinContent( 4 );
	err_250   = h_sigRegion->GetBinError( 1 );
	err_300   = h_sigRegion->GetBinError( 2 );
	err_350   = h_sigRegion->GetBinError( 3 );
	err_400   = h_sigRegion->GetBinError( 4 );

	// Print LaTeX table line
	printf("%28s  &   %8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f   &   ", printNames.at(i).Data(), yield_250, err_250, yield_300, err_300);
	printf("%8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f  \\\\\n", yield_350, err_350, yield_400, err_400);
	if( sampleNames.at(i)=="stop425" ) cout << "\\hline" << endl;
  } // end loop over backgrounds
  cout << "\\hline" << endl;

  // Get total yields
  yield_250 = h_totals_sregion->GetBinContent( 1 );
  yield_300 = h_totals_sregion->GetBinContent( 2 );
  yield_350 = h_totals_sregion->GetBinContent( 3 );
  yield_400 = h_totals_sregion->GetBinContent( 4 );
  err_250   = h_totals_sregion->GetBinError( 1 );
  err_300   = h_totals_sregion->GetBinError( 2 );
  err_350   = h_totals_sregion->GetBinError( 3 );
  err_400   = h_totals_sregion->GetBinError( 4 );

  // Print LaTeX table line
  printf("%28s  &   %8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f   &   ", "Total background", yield_250, err_250, yield_300, err_300);
  printf("%8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f  \\\\\n", yield_350, err_350, yield_400, err_400);

  cout << "\\hline \\hline" << endl;


  /////////////////////////////////////////////////////////////////////////////////
  // Make the lower half of the table, with the rows = final states

  for( uint i=0; i<bkgNames.size()-1; i++ ) {
	// Grab the yields and errors by background type
	yield_250 = h_bgtype[0]->GetBinContent( i+1 );
	yield_300 = h_bgtype[1]->GetBinContent( i+1 );
	yield_350 = h_bgtype[2]->GetBinContent( i+1 );
	yield_400 = h_bgtype[3]->GetBinContent( i+1 );
	err_250   = h_bgtype[0]->GetBinError( i+1 );
	err_300   = h_bgtype[1]->GetBinError( i+1 );
	err_350   = h_bgtype[2]->GetBinError( i+1 );
	err_400   = h_bgtype[3]->GetBinError( i+1 );

	// Print LaTeX table line
	printf( "%28s   &   %8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f   &   ", bkgNames.at(i).Data(), yield_250, err_250, yield_300, err_300 );
	printf( "%8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f   \\\\\n", yield_350, err_350, yield_400, err_400 );
  }

  cout << "\\end{tabular}\n\n" << endl;



  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now do the second table, for the high Delta-M region

  // Make the upper part of the table, with rows = different samples, columns = signal regions

  double yield_500, err_500;

  cout << "\\begin{tabular}{ l | c c c c c }" << endl;
  cout << "High $\\Delta M$   &  250 $\\leq E_T^{miss} <$ 300  &   300 $\\leq E_T^{miss} <$ 350  &   350 $\\leq E_T^{miss} <$ 400  &   400 $\\leq E_T^{miss} <$ 500 & $E_T^{miss} \\geq$ 500  \\\\ \\hline" << endl;

  for( uint i=0; i<sampleNames.size(); i++ ) {
	// Read in yields and errors from stored histograms
	TH1D* h_sigRegion = (TH1D*)infile->Get("sregion_" + sampleNames.at(i));
	// h_totals_sregion->Add(h_sigRegion); // already summed up backgrounds in the first table; don't need to repeat

	yield_250 = h_sigRegion->GetBinContent( 5 );
	yield_300 = h_sigRegion->GetBinContent( 6 );
	yield_350 = h_sigRegion->GetBinContent( 7 );
	yield_400 = h_sigRegion->GetBinContent( 8 );
	yield_500 = h_sigRegion->GetBinContent( 9 );
	err_250   = h_sigRegion->GetBinError( 5 );
	err_300   = h_sigRegion->GetBinError( 6 );
	err_350   = h_sigRegion->GetBinError( 7 );
	err_400   = h_sigRegion->GetBinError( 8 );
	err_500   = h_sigRegion->GetBinError( 9 );

	// Print LaTeX table line
	printf("%28s  &   %8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f   &   ", printNames.at(i).Data(), yield_250, err_250, yield_300, err_300);
	printf("%8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f &   %8.3f $\\pm$ %6.3f  \\\\\n", yield_350, err_350, yield_400, err_400, yield_500, err_500);
	if( sampleNames.at(i)=="stop425" ) cout << "\\hline" << endl;
  } // end loop over backgrounds
  cout << "\\hline" << endl;

  // Get total yields
  yield_250 = h_totals_sregion->GetBinContent( 5 );
  yield_300 = h_totals_sregion->GetBinContent( 6 );
  yield_350 = h_totals_sregion->GetBinContent( 7 );
  yield_400 = h_totals_sregion->GetBinContent( 8 );
  yield_500 = h_totals_sregion->GetBinContent( 9 );
  err_250   = h_totals_sregion->GetBinError( 5 );
  err_300   = h_totals_sregion->GetBinError( 6 );
  err_350   = h_totals_sregion->GetBinError( 7 );
  err_400   = h_totals_sregion->GetBinError( 8 );
  err_500   = h_totals_sregion->GetBinError( 9 );

  // Print LaTeX table line
  printf("%28s  &   %8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f   &   ", "Total background", yield_250, err_250, yield_300, err_300);
  printf("%8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f &   %8.3f $\\pm$ %6.3f  \\\\\n", yield_350, err_350, yield_400, err_400, yield_500, err_500);

  cout << "\\hline \\hline" << endl;


  /////////////////////////////////////////////////////////////////////////////////
  // Make the lower half of the table, with the rows = final states

  for( uint i=0; i<bkgNames.size()-1; i++ ) {
	// Grab the yields and errors by background type
	yield_250 = h_bgtype[4]->GetBinContent( i+1 );
	yield_300 = h_bgtype[5]->GetBinContent( i+1 );
	yield_350 = h_bgtype[6]->GetBinContent( i+1 );
	yield_400 = h_bgtype[7]->GetBinContent( i+1 );
	yield_500 = h_bgtype[8]->GetBinContent( i+1 );
	err_250   = h_bgtype[4]->GetBinError( i+1 );
	err_300   = h_bgtype[5]->GetBinError( i+1 );
	err_350   = h_bgtype[6]->GetBinError( i+1 );
	err_400   = h_bgtype[7]->GetBinError( i+1 );
	err_500   = h_bgtype[8]->GetBinError( i+1 );

	// Print LaTeX table line
	printf( "%28s   &   %8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f   &   ", bkgNames.at(i).Data(), yield_250, err_250, yield_300, err_300 );
	printf( "%8.3f $\\pm$ %6.3f   &   %8.3f $\\pm$ %6.3f &   %8.3f $\\pm$ %6.3f   \\\\\n", yield_350, err_350, yield_400, err_400, yield_500, err_500 );
  }

  cout << "\\end{tabular}\n\n" << endl;


}
