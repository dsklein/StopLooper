#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"

#include "analysis.h"
#include "sample.h"

using namespace std;

double calculateMaxdiff( int bin, TFile* histfile, const double& nominal, vector<TString> varNames );
void printSystTable( vector<TString> sigRegions, map<TString,vector<double> > uncertainties );


void makeDataCards( analysis* srAnalysis, analysis* lostlepAnalysis = NULL, analysis* onelepwAnalysis = NULL ) {

	// Do the basic setup stuff

	vector<TString> bkgs = { "zNuNu", "dilep", "top1l", "W1l" }; // Eventually pull this from the analysis object
	const int nBkgs = bkgs.size();

	if( srAnalysis->GetNsignals() < 1 ) {
		cout << "\nError in makeDataCards.cc: Need at least one signal sample!" << endl;
		return;
	}

	vector<TString> samples = bkgs;
	samples.insert( samples.begin(), "signal" );
	const int nSamples = samples.size();

	vector<TString> sigRegions = srAnalysis->GetSigRegionLabelsAll();
	const int nSigRegs = sigRegions.size();

	int nVars_ll, nVars_1lw;
	map<TString,vector<TString> > systMap_sr, systMap_ll, systMap_1lw;

	systMap_sr = srAnalysis->GetSystMap();
	if( lostlepAnalysis ) {
		nVars_ll  = lostlepAnalysis->GetSystematics(true).size();
		systMap_ll = lostlepAnalysis->GetSystMap();
	}
	if( onelepwAnalysis ) {
		nVars_1lw   = onelepwAnalysis->GetSystematics(true).size();
		systMap_1lw = onelepwAnalysis->GetSystMap();
	}

	map<TString,vector<double> > lostlep_uncert, onelepw_uncert;

	// Open files containing background yields and uncertainties
	TFile* yieldFile   = new TFile( srAnalysis->GetPlotFileName(), "READ" );
	TFile* lostlepFile = new TFile( "lostlepEstimates.root", "READ" );
	TFile* onelepwFile = new TFile( "onelepwEstimates.root", "READ" );
	for( TFile* thisFile : {yieldFile, lostlepFile, onelepwFile} ) {
		if( thisFile->IsZombie() ) {
			cout << "Error in makeDataCards! Couldn't open file" << thisFile->GetName() << "!" << endl;
			return;
		}
	}

	TH1D* h_lostLep = (TH1D*)lostlepFile->Get("lostLepBkg");
	TH1D* h_onelepw = (TH1D*)onelepwFile->Get("onelepwBkg");

	//////////////////////////////////////////////////////////////////////////////
	// Loop over signal regions, making a datacard for each SR and each mass point
	//////////////////////////////////////////////////////////////////////////////

	for( int reg=1; reg<=nSigRegs; reg++ ) {

		cout << "Writing data cards for signal region " << sigRegions.at(reg-1) << endl;

		TH1D* h_bkgYield  = (TH1D*)yieldFile->Get( "evttype_"+sigRegions.at(reg-1) );
		TH2D* h_sigYield  = (TH2D*)yieldFile->Get( "sigyields_"+sigRegions.at(reg-1) );
		TH2D* h_sigContam = (TH2D*)lostlepFile->Get( "sigContam_"+sigRegions.at(reg-1) );

		// Subtract signal contamination from signal MC prediction.
		// Had a LONG discussion with John and FKW about why they want us to do it this way
		h_sigYield->Add( h_sigContam, -1. );

		// Get bin sizes for the signal yield histogram
		double binwidthX = h_sigYield->GetXaxis()->GetBinWidth(1);
		double binwidthY = h_sigYield->GetYaxis()->GetBinWidth(1);


		///////////////////////////////////////////////////////////////////////////////////////
		// Begin constructing a template datacard for this particular signal region.

		vector<TString> cardlines, uncertlines;
		TString tmpstr;


		// Top matter

		cardlines.push_back( Form( "# Data card for signal region %d (%s)\n", reg, sigRegions.at(reg-1).Data() ) );
		cardlines.push_back( "### Placeholder line. SUSY masses will go here." );
		cardlines.push_back( "---\n" );
		cardlines.push_back( "imax 1  number of channels\n" );

		tmpstr = Form( "jmax %d  number of backgrounds (", nBkgs );
		for (TString bkgName : bkgs ) tmpstr += Form( "%s, ",  bkgName.Data() );   // Will usually be 4 (znunu, 2l, 1ltop, 1lw)
		tmpstr += ")\n";
		cardlines.push_back( tmpstr );
		cardlines.push_back( "### Placeholder line. Number of uncertainties will go here.\n" );
		cardlines.push_back( "---\n" );

		cardlines.push_back( "# Now list the number of events observed (or zero if no data)\n" );
		cardlines.push_back( Form("bin %d\n", reg) );
		cardlines.push_back( Form( "observation %d\n", int(h_bkgYield->GetBinContent(1)) ) );
		cardlines.push_back( "---\n" );

		cardlines.push_back(  "# Now list the expected events (i.e. Monte Carlo yield) for signal and all backgrounds in our particular bin\n" );
		cardlines.push_back(  "# Second 'process' line should be zero for signal, and a positive integer for each background\n" );

		tmpstr = "bin      ";
		for( TString sample : samples ) tmpstr += Form( " %10i", reg );  // Write the signal region number once for each sample (signal and bkg)
		tmpstr += "\n";
		cardlines.push_back( tmpstr );

		tmpstr = "process  ";
		for( TString sample : samples ) tmpstr += Form( " %10s", sample.Data() );  // Write the name of each sample (sig & bkg)
		tmpstr += "\n";
		cardlines.push_back( tmpstr );

		tmpstr = "process  ";
		for( int j=0; j<nSamples; j++ ) tmpstr += Form( " %10i", j ); // Write a number for each sample (0=signal, positive integers for bkgs)
		tmpstr += "\n";
		cardlines.push_back( tmpstr );


		// Precalculate the yields for each background process

		double yield;
		TString bkgRates;

		for( int i=1; i<nSamples; i++ ) {
			if(      i==2 ) yield = h_lostLep->GetBinContent(reg); // Pull 2l yield from lostLepton estimate histogram
			else if( i==4 ) yield = h_onelepw->GetBinContent(reg); // Pull 1lw yield from 1l-from-w estimate histogram
			else            yield = h_bkgYield->GetBinContent(i+2);
			bkgRates += Form( " %10f", yield );
		}
		bkgRates += "\n";
		cardlines.push_back( "### Placeholder line. Observed rate of background processes will go here." );

		cardlines.push_back( "---\n" );
		cardlines.push_back( "# Now we list the independent sources of uncertainty (syst. and stat. error), and which samples they affect\n" );
		cardlines.push_back( "---\n" );


		/////////////////////////////////////////////////////////////////////////
		// Make a similar template for the rows that store all the uncertainties

		// Generate a row for the statistical uncertainty on each sample
		//  (except signal, which is handled later per-mass-point, and 1l-from-top, which is covered by the 100% uncertainty)
		uncertlines.push_back( "### Placeholder line. Stat uncertainty on signal yield will go here." );
		for( int sampleIdx=1; sampleIdx<nSamples; sampleIdx++ ) {

			if( sampleIdx == 3 ) continue;

			char statname[25];
			sprintf( statname, "Stat%s%d", samples.at(sampleIdx).Data(), reg );
			tmpstr = Form( "%-18s  lnN ", statname );

			double statErr = 1.0 + (h_bkgYield->GetBinError(sampleIdx+2) / h_bkgYield->GetBinContent(sampleIdx+2) );
			if(      sampleIdx==2 ) statErr = 1.0 + ( h_lostLep->GetBinError(reg) / h_lostLep->GetBinContent(reg) );  // Pull 2l stat error from lostLepton estimate histogram
			else if( sampleIdx==4 ) statErr = 1.0 + ( h_onelepw->GetBinError(reg) / h_onelepw->GetBinContent(reg) );  // Pull 1lw stat error from 1l-from-w estimate histogram

			if( std::isnan(statErr) ) statErr = 1.0; // Protection against nan

			for( int j=0; j<nSamples; j++ ) {
				if( j == sampleIdx )  tmpstr += Form( "  %8.6f  ", statErr);
				else tmpstr += "     -      ";
			}
			tmpstr += "\n";
			uncertlines.push_back( tmpstr );
		}


		// Write out a dummy systematic uncertainty for each of the backgrounds that doesn't have an actual systematic calculation
		//   For now, skip the signal sample (don't give it a systematic uncertainty)
		for( int sampleIdx = 1; sampleIdx<nSamples; sampleIdx++ ) {

			double systErr = 0.;
			char systname[25];

			if(      sampleIdx == 2 && nVars_ll  > 0 ) continue; // Don't use dummy systematic for ll background if we have actual systematics
			else if( sampleIdx == 4 && nVars_1lw > 0 ) continue; // Same with 1l-from-W background
			else if( sampleIdx == 3 ) {
				systErr = 2.0; // 100% systematic on 1l from top
				sprintf( systname, "Flat%s%d", samples.at(sampleIdx).Data(), reg );
			}
			else {
				systErr = 1.3; // Flat 30% systematic for now
				sprintf( systname, "Flat%s", samples.at(sampleIdx).Data() );
			}

			tmpstr = Form( "%-18s  lnN ", systname );

			for( int j=0; j<nSamples; j++ ) {
				if( j == sampleIdx ) tmpstr += Form( "  %8.6f  ", systErr);
				else tmpstr +=  "     -      " ;
			}
			tmpstr += "\n";
			uncertlines.push_back( tmpstr );
		}

		////////////////////////////////////////////////////////////////
		// Make the systematics section of the datacard

		// Start with an empty systematics table...
		map<TString,vector<TString> > syst_holder;
		vector<TString> emptySystLine = {"     -      ", "     -      ", "     -      ", "     -      ", "     -      " };
		for( auto& iter : systMap_ll )  syst_holder[iter.first] = emptySystLine;
		for( auto& iter : systMap_1lw ) syst_holder[iter.first] = emptySystLine;

		// Calculate the lost lepton systematics, and populate the datacard and the printable systematic tables
		if( nVars_ll > 0 ) {
			double nominal = h_lostLep->GetBinContent(reg);
			if( nominal < 0.0000000001 ) nominal = 0.0000000001;  // Protection for when yields are zero or negative
			for( auto& iter : systMap_ll ) {                      // Loop over all systematics for the ll background
				double maxdiff = calculateMaxdiff( reg, lostlepFile, nominal, iter.second ); // Calculate the biggest variation
				syst_holder[iter.first].at(2) = Form( "  %8.6f  ", 1.0 + maxdiff/nominal );  // Populate the appropriate cell in the datacard systematics table
				lostlep_uncert[iter.first].push_back( maxdiff/nominal );     // and also add that number to the printable lost lepton systematics table
			}
		}

		// Same with 1l-from-W systematics
		if( nVars_1lw > 0 ) {
			double nominal = h_onelepw->GetBinContent(reg);
			if( nominal < 0.0000000001 ) nominal = 0.0000000001;
			for( auto& iter : systMap_1lw ) {
				double maxdiff = calculateMaxdiff( reg, onelepwFile, nominal, iter.second );
				syst_holder[iter.first].at(4) = Form( "  %8.6f  ", 1.0 + maxdiff/nominal );
				onelepw_uncert[iter.first].push_back( maxdiff/nominal );
			}
		}

		// Now assemble the systematics table, and add it to our datacard template
		for( auto& iter : syst_holder ) {
			tmpstr = Form( "Syst%-14s  lnN ", iter.first.Data() );
			for( TString cell : iter.second ) tmpstr += cell;
			tmpstr += "\n";
			uncertlines.push_back( tmpstr );
		}

		// Count number of uncertainties, and insert appropriate row into datacard template
		int nUncerts = uncertlines.size();
		cardlines.at(5) = Form( "kmax %d  number of uncertainties\n", nUncerts );



		///////////////////////////////////////////////////////////////////////////////////////////
		// Loop over signal mass points, to get the information that's specific to each mass point

		for( int xbin=1; xbin<=h_sigYield->GetNbinsX(); xbin++ ) {
			for( int ybin=1; ybin<=h_sigYield->GetNbinsY(); ybin++ ) {

				// Round bin centers to sensible numbers (nearest integer multiple of the bin width)
				int stopmass = binwidthX * round( h_sigYield->GetXaxis()->GetBinCenter(xbin) / binwidthX );
				int lspmass  = binwidthY * round( h_sigYield->GetYaxis()->GetBinCenter(ybin) / binwidthY );

				double sigYield = h_sigYield->GetBinContent( xbin, ybin );
				double sigError = h_sigYield->GetBinError(   xbin, ybin );
				if( fabs(sigYield) < 0.000001 ) continue;
				if( sigYield < 0. ) {
					// cout << "Warning in makeDataCards: Mass point (" << stopmass << "," << lspmass
					// 	   << ") has negative signal yield after subtracting contamination. Setting yield to zero." << endl;
					sigYield = 0.000000001;
					sigError = 0.;
				}

				// Generate the second line of the datacard, with the susy masses
				cardlines.at(1) = Form( "# Stop mass = %d, LSP mass = %d\n", stopmass, lspmass );

				// Generate line 17, with the "rate" for signal and backgrounds
				tmpstr = "rate     ";
				tmpstr += Form( " %10f", sigYield );
				tmpstr += bkgRates;
				cardlines.at(16) = tmpstr;

				// Generate the line that holds the stat uncertainty on signal
				char statname[25];
				sprintf( statname, "Stat%s%d", samples.at(0).Data(), reg );
				tmpstr = Form( "%-18s  lnN ", statname );

				double statErr = 1.0 + ( sigError / sigYield );
				tmpstr += Form( "  %8.6f  ", statErr);
				for( int j=1; j<nSamples; j++ ) tmpstr += "     -      ";
				tmpstr += "\n";
				uncertlines.at(0) = tmpstr;


				///////////////////////////////////////////////////////////////////
				// Now let's actually make a datacard!

				TString fileName = Form( "datacards/datacard_%s_T2tt_%d_%d.txt", sigRegions.at(reg-1).Data(), stopmass, lspmass );

				// Open file
				FILE * outfile;
				outfile = fopen( fileName.Data(), "w" );

				// Write out each line
				for( TString line : cardlines   ) fputs( line.Data(), outfile );
				for( TString line : uncertlines ) fputs( line.Data(), outfile );

				// Close file
				fprintf( outfile,  "---\n" );
				fclose(outfile);

			} // End loop over y bins (LSP masses)
		} // End loop over x bins (stop masses)

	} // End loop over signal regions



	// Print table of systematics for the dilepton background estimate and the 1l-from-W estimate
	if( nVars_ll > 0 ) {
		printf( "\n\nSystematics on dilepton background estimate\n\n" );
		printSystTable( sigRegions, lostlep_uncert );
	}
	if( nVars_1lw > 0 ) {
		printf( "\n\nSystematics on 1l-from-W background estimate\n\n" );
		printSystTable( sigRegions, onelepw_uncert );
	}


	lostlepFile->Close();
	onelepwFile->Close();
	yieldFile->Close();
	delete lostlepFile;
	delete onelepwFile;
	delete yieldFile;
}




double calculateMaxdiff( int bin, TFile* histfile, const double& nominal, vector<TString> varNames ) {

	if( nominal < 0.00000001 ) return 0.;
	double maxdiff = 0.;

	for( TString varName : varNames ) {
		TH1D* h_tmp = (TH1D*)histfile->Get( "variation_" + varName );
		if( h_tmp == 0 ) {
			cout << "Warning in makeDataCards: Couldn't find histogram variation_" << varName << " in file " << histfile->GetName() << "!" << endl;
			continue;
		}
		maxdiff = max( maxdiff, fabs(nominal - h_tmp->GetBinContent(bin)) );
	}

	return maxdiff;
}


void printSystTable( vector<TString> sigRegions, map<TString,vector<double> > uncertainties ) {

	printf( "\n\\begin{tabular}{ | l |" );
	for( TString regName : sigRegions ) printf( " c |" );
	printf( " }\n" );
	printf( "\\hline\n" );

	printf( "Systematic " );
	for( TString regName : sigRegions ) printf( "& %s ", regName.Data() );
	printf( " \\\\ \\hline\n" );

	for( auto& iter : uncertainties ) {
		printf( "%10s ", iter.first.Data() );
		for( double uncert : iter.second ) printf( "& %4.1f\\%% ", uncert*100. );
		printf( " \\\\\n" );
	}

	printf( "\\hline\n" );
	printf( "\\end{tabular}\n\n" );

}
