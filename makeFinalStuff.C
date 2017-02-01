#include <iostream>

#include "TColor.h"
#include "TFile.h"
#include "TH1F.h"

// #include "analysis.h"
// #include "sample.h"
// #include "sigRegion.h"

#include "dataMCplotMaker.h"

using namespace std;



TH1F* shortenHisto( TH1F* hist );
TH1F* makeErrorHisto( TH1F* hist, double e1, double e2, double e3, double e4 );



void makeFinalStuff() {

	// Open files
	TFile* f_sr    = new TFile( "plots.root", "READ" );
	TFile* f_cr2l  = new TFile( "lostlepEstimates.root", "READ" );
	TFile* f_cr0b  = new TFile( "onelepwEstimates.root", "READ" );
	TFile* f_znunu = new TFile( "ZNuNu_BkgEst.root", "READ" );
	TFile* f_1ltop = new TFile( "oneleptopEstimates.root", "READ" );

	// Get histograms
	TH1F* h_data    = (TH1F*)f_sr->Get("srYields_data");
	TH1F* h_lostlep = (TH1F*)f_cr2l->Get("lostLepBkg");
	TH1F* h_onelepw = (TH1F*)f_cr0b->Get("onelepwBkg");
	TH1F* h_znunu_raw = (TH1F*)f_znunu->Get("yield");
	TH1F* h_onelept = (TH1F*)f_1ltop->Get("srMC");

	// Correct the binning on the Znunu histogram
	TH1F* h_znunu = (TH1F*)h_lostlep->Clone("znunu");
	h_znunu->Reset();
	for( int i=1; i<=31; i++ ) {
		h_znunu->SetBinContent(i, h_znunu_raw->GetBinContent(i));
		h_znunu->SetBinError(i, h_znunu_raw->GetBinError(i));
	}

	// Set 100% uncertainty on 1l from top for the purposes of the table (will transfer this uncertainty to systematic later)
	for( int i=1; i<=31; i++ ) h_onelept->SetBinError( i, h_onelept->GetBinContent(i) );

	// Scale to the appropriate luminosity
	// h_lostlep->Scale( 18.1 / 36.46 );
	// h_onelepw->Scale( 18.1 / 36.46 );
	// h_znunu->Scale(   18.1 / 36.46 );
	// h_onelept->Scale( 18.1 / 36.46 );

	// Make histo of total bkg yield
	TH1F* h_totalbkg = (TH1F*)h_lostlep->Clone("total_bkg");
	h_totalbkg->Add( h_onelepw );
	h_totalbkg->Add( h_znunu );
	h_totalbkg->Add( h_onelept );

	vector<TH1F*> histos = {h_lostlep, h_onelept, h_onelepw, h_znunu, h_totalbkg, h_data};

	// Make yield table
	printf( "\\begin{tabular}{| l | c | c | c | c | c | c | }\n" );
	printf( "\\hline\n" );
	printf( "Region  &    Lost lepton  &  1$\\ell$ (top)  &  1$\\ell$ (W)  &  $Z\\rightarrow\\nu\\nu$  & Total background  & Data \\\\\n" );
	printf( "\\hline\n" );
	for( int i=28; i<=31; i++ ) {
		printf( "%12s  ", h_data->GetXaxis()->GetBinLabel(i) );
		for( TH1F* hist : histos ) printf( " &  %5.2f $\\pm$ %5.2f ", hist->GetBinContent(i), hist->GetBinError(i) );
		printf( " \\\\\n" );
	}
	printf( "\\hline\n" );
	printf( "\\end{tabular}\n" );




	//////////////////////////////////////
	// Now make plots

	// Move 100% uncertainty on 1l top from stat to systematic
	for( int i=1; i<=31; i++ ) h_onelept->SetBinError( i, 0. );

	TH1F* h_data_shortened = shortenHisto( h_data );
	TH1F* h_lostlep_short = shortenHisto( h_lostlep );
	TH1F* h_onelept_short = shortenHisto( h_onelept );
	TH1F* h_onelepw_short = shortenHisto( h_onelepw );
	TH1F* h_znunu_short   = shortenHisto( h_znunu );


	// vector<TH1F*> bkgs_big = {h_lostlep, h_onelept, h_onelepw, h_znunu};
	vector< pair<TH1F*,TH1F*> > bkgs;
	vector<string> bkg_titles = {"Lost lepton", "1l (top)", "1l (W)", "Z#rightarrow#nu#nu"};
	vector<short> colors = {kCyan-3, kRed-7, kOrange-2, kMagenta-5};
	vector<TH1F*> sigs;
	vector<string> sig_titles;


	// Manually input systematic errors for each bkg component
	TH1F* h_error_ll = makeErrorHisto( h_lostlep_short, 0.141, 0.236, 0.538, 0.519 );
	TH1F* h_error_1lt = makeErrorHisto( h_onelept_short, 1., 1., 1., 1. );
	TH1F* h_error_1lw = makeErrorHisto( h_onelepw_short, 0.419, 0.404, 0.658, 0.638 );
	TH1F* h_error_znunu = makeErrorHisto( h_znunu_short, 0.1454, 0.1785, 0.2015, 0.2430 );

	// Make a vector of pairs of histograms and errors
	bkgs.push_back( make_pair(h_lostlep_short, h_error_ll) );
	bkgs.push_back( make_pair(h_onelept_short, h_error_1lt) );
	bkgs.push_back( make_pair(h_onelepw_short, h_error_1lw) );
	bkgs.push_back( make_pair(h_znunu_short, h_error_znunu) );






	dataMCplotMaker( h_data_shortened,
	                 bkgs,
	                 bkg_titles,
	                 "",
	                 "",
	                 "--energy 13 --lumi 36.46 --type Preliminary --outOfFrame --xAxisLabel Region --yAxisLabel Events --noXaxisUnit --outputName finalPlot.pdf",
	                 sigs,
	                 sig_titles,
	                 colors );

	dataMCplotMaker( h_data_shortened,
	                 bkgs,
	                 bkg_titles,
	                 "",
	                 "",
	                 "--energy 13 --lumi 18.1 --type Preliminary --outOfFrame --xAxisLabel Region --yAxisLabel Events --noXaxisUnit --outputName finalPlot.eps",
	                 sigs,
	                 sig_titles,
	                 colors );

}






TH1F* shortenHisto( TH1F* hist ) {

	char newname[30];
	sprintf( newname, "%s_copy", hist->GetName() );

	TH1F* h_out = new TH1F( newname, hist->GetTitle(), 4, 0.5, 4.5 );
	int nbins = hist->GetNbinsX();
	int bindiff = nbins - 4;

	for( int i=1; i<=4; i++ ) {
		h_out->SetBinContent( i, hist->GetBinContent(bindiff+i) );
		h_out->SetBinError( i, hist->GetBinError(bindiff+i) );
	}

	TAxis* ax = h_out->GetXaxis();
	ax->SetBinLabel( 1, "Corridor 250-350" );
	ax->SetBinLabel( 2, "Corridor 350-450" );
	ax->SetBinLabel( 3, "Corridor 450-550" );
	ax->SetBinLabel( 4, "Corridor 550+" );

	return h_out;
}

TH1F* makeErrorHisto( TH1F* hist, double e1, double e2, double e3, double e4 ) {

	char newname[30];
	sprintf( newname, "%s_error", hist->GetName() );

	TH1F* h_out = (TH1F*)hist->Clone(newname);
	h_out->SetBinError( 1, e1*hist->GetBinContent(1) );
	h_out->SetBinError( 2, e2*hist->GetBinContent(2) );
	h_out->SetBinError( 3, e3*hist->GetBinContent(3) );
	h_out->SetBinError( 4, e4*hist->GetBinContent(4) );
	// h_out->Multiply( hist );

	return h_out;
}
