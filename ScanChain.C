// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>
#include <cmath>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TVector3.h"
#include "Math/VectorUtil.h"

// CMS3
#include "CMS3.cc"

using namespace std;
using namespace tas;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

double myLumi = 5.;

int ScanChain( TChain* chain, string sampleName = "default", int nEvents = -1, bool fast = true) {

  cout << "\nSample: " << sampleName << endl;

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  /////////////////////////////////////////////////////////
  // Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  TH1::SetDefaultSumw2();


  TH1D* h_bgtype[9];

  TH1D *h_mt[9];
  TH1D *h_met[9];
  TH1D *h_mt2w[9];
  TH1D *h_chi2[9];
  TH1D *h_htratio[9];
  TH1D *h_mindphi[9];
  TH1D *h_ptb1[9];
  TH1D *h_drlb1[9];
  TH1D *h_ptlep[9];
  TH1D *h_metht[9];
  TH1D *h_dphilw[9];
  TH1D *h_njets[9];
  TH1D *h_nbtags[9];

  double met_min[9]  = {250., 300., 350., 400.,   250., 300., 350., 400., 500.};
  double met_max[9]  = {300., 350., 400., 99999., 300., 350., 400., 500., 99999.};
  double mt2w_min[9] = {  0.,   0.,   0.,   0.,   200., 200., 200., 200., 200.};
  double mt2w_max[9] = {200., 200., 200., 200., 99999., 99999., 99999., 99999., 99999.};
  string regNames[9] = {"low250", "low300", "low350", "low400", "high250", "high300", "high350", "high400", "high500"};


  for( int i=0; i<9; i++ ) {

	h_bgtype[i]   = new TH1D( Form( "bkgtype_%s_%s" , sampleName.c_str(), regNames[i].c_str()), "Yield by background type", 4, 0.5, 4.5);
	h_mt[i]       = new TH1D( Form( "mt_%s_%s"      , sampleName.c_str(), regNames[i].c_str()),	"Transverse mass",			80, 0, 800);
	h_met[i]      = new TH1D( Form( "met_%s_%s"     , sampleName.c_str(), regNames[i].c_str()),	"MET",						50, 0, 1000);
	h_mt2w[i]     = new TH1D( Form( "mt2w_%s_%s"    , sampleName.c_str(), regNames[i].c_str()),	"MT2W",						50, 0, 500);
	h_chi2[i]     = new TH1D( Form( "chi2_%s_%s"    , sampleName.c_str(), regNames[i].c_str()),	"Hadronic #chi^{2}",		50, 0, 15);
	h_htratio[i]  = new TH1D( Form( "htratio_%s_%s" , sampleName.c_str(), regNames[i].c_str()),	"H_{T} ratio",				50, 0, 1);
	h_mindphi[i]  = new TH1D( Form( "mindphi_%s_%s" , sampleName.c_str(), regNames[i].c_str()),	"min #Delta#phi(j12,MET)",	50, 0, 4);
	h_ptb1[i]     = new TH1D( Form( "ptb1_%s_%s"    , sampleName.c_str(), regNames[i].c_str()),	"p_{T} (b1)",				100, 0, 500);
	h_drlb1[i]    = new TH1D( Form( "drlb1_%s_%s"   , sampleName.c_str(), regNames[i].c_str()),	"#DeltaR (lep, b1)",		50, 0, 5);
	h_ptlep[i]    = new TH1D( Form( "ptlep_%s_%s"   , sampleName.c_str(), regNames[i].c_str()),	"p_{T} (lep)",				100, 0, 500);
	h_metht[i]    = new TH1D( Form( "metht_%s_%s"   , sampleName.c_str(), regNames[i].c_str()),	"MET/sqrt(HT)",				50, 0, 100);
	h_dphilw[i]   = new TH1D( Form( "dphilw_%s_%s"  , sampleName.c_str(), regNames[i].c_str()),	"#Delta#phi (lep,W)",		50, 0, 3.5);
	h_njets[i]    = new TH1D( Form( "njets_%s_%s"   , sampleName.c_str(), regNames[i].c_str()), "Number of jets",           16, -0.5, 15.5);
	h_nbtags[i]   = new TH1D( Form( "nbtags_%s_%s"  , sampleName.c_str(), regNames[i].c_str()), "Number of b-tags",         7, -0.5, 6.5);

	h_bgtype[i]->SetDirectory(rootdir);

	h_mt[i]->SetDirectory(rootdir);
	h_met[i]->SetDirectory(rootdir);
	h_mt2w[i]->SetDirectory(rootdir);
	h_chi2[i]->SetDirectory(rootdir);
	h_htratio[i]->SetDirectory(rootdir);
	h_mindphi[i]->SetDirectory(rootdir);
	h_ptb1[i]->SetDirectory(rootdir);
	h_drlb1[i]->SetDirectory(rootdir);
	h_ptlep[i]->SetDirectory(rootdir);
	h_metht[i]->SetDirectory(rootdir);
	h_dphilw[i]->SetDirectory(rootdir);
	h_njets[i]->SetDirectory(rootdir);
	h_nbtags[i]->SetDirectory(rootdir);

	TAxis* axis = h_bgtype[i]->GetXaxis();
	axis->SetBinLabel( 1, "1lep" );
	axis->SetBinLabel( 2, "2+lep" );
	axis->SetBinLabel( 3, "ZtoNuNu" );
	axis->SetBinLabel( 4, "Other" );

  }

  TH1D *h_sigRegion = new TH1D( Form("sregion_%s", sampleName.c_str()), "Yield by signal region", 9, 0.5, 9.5);
  h_sigRegion->SetDirectory(rootdir);


  float yield_total = 0;
  float yield_vtx = 0;
  float yield_1goodlep = 0;
  float yield_2lepveto = 0;
  float yield_lepSel = 0;
  float yield_trkVeto = 0;
  float yield_tauVeto = 0;
  float yield_4jets = 0;
  float yield_1bjet = 0;
  float yield_METcut = 0;
  float yield_MTcut = 0;
  float yield_dPhi = 0;
  float yield_chi2 = 0;

  int yGen_total = 0;
  int yGen_vtx = 0;
  int yGen_1goodlep = 0;
  int yGen_2lepveto = 0;
  int yGen_lepSel = 0;
  int yGen_trkVeto = 0;
  int yGen_tauVeto = 0;
  int yGen_4jets = 0;
  int yGen_1bjet = 0;
  int yGen_METcut = 0;
  int yGen_MTcut = 0;
  int yGen_dPhi = 0;
  int yGen_chi2 = 0;


  /////////////////////////////////////////////////////////////////////

  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("t");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    cms3.Init(tree);
    
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      cms3.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      CMS3::progress( nEventsTotal, nEventsChain );

	  ////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Analysis Code
	  // ---------------------------------------------------------------------------------------------------//


	  ///////////////////////////////////////////////////////////////
	  // Special filters to more finely categorize background events
	  if(      sampleName == "Wb"    && genbs_p4().size() < 1 ) continue;  //Make sure there are b-quarks in the "W+b" events
	  else if( sampleName == "Wucsd" && genbs_p4().size() > 0 ) continue;  //Make sure there are no b-quarks in the "W+ucsd" events
	  else if( sampleName == "tt2l"  && genlepsfromtop() != 2 ) continue;  //Require 2 leps from top in "tt2l" events
	  else if( sampleName == "tt1l"  && genlepsfromtop() != 1 ) continue;  //Require 1 lep from top in "tt1l" events


	  // Count the number of events processed
	  yield_total += myLumi*scale1fb();
	  yGen_total++;

	  // First vertex must be good
	  if( firstGoodVtxIdx() != 0 ) continue;
	  yield_vtx += myLumi*scale1fb();
	  yGen_vtx++;

	  // Must have exactly 1 good lepton
	  if( ngoodleps() != 1 ) continue;
	  yield_1goodlep += myLumi*scale1fb();
	  yGen_1goodlep++;

	  // Second lepton veto
	  if( nvetoleps() > 1 && ROOT::Math::VectorUtil::DeltaR( lep1_p4(), lep2_p4() ) > 0.01 ) continue;
	  yield_2lepveto += myLumi*scale1fb();
	  yGen_2lepveto++;

	  // Must pass lepton selections
	  if( lep1_is_el() ) {
		if( lep1_pt() <= 40. ) continue;
		if( fabs(lep1_eta()) >= 2.1 ) continue;
		if( !lep1_is_phys14_medium_noIso() ) continue;
		if( lep1_miniRelIsoDB() >= 0.1 ) continue;
	  }
	  else if( lep1_is_mu() ) {
		if( lep1_pt() <= 30. ) continue;
		if( fabs(lep1_eta()) >= 2.1 ) continue;
		if( !lep1_is_muoid_medium() ) continue;
		if( fabs(lep1_d0()) >= 0.02 ) continue;
		if( fabs(lep1_dz()) >= 0.1  ) continue;
		if( lep1_miniRelIsoDB() >= 0.1 ) continue;
	  }
	  yield_lepSel += myLumi*scale1fb();
	  yGen_lepSel++;

	  // Track veto
	  if( !PassTrackVeto_v3() ) continue;
	  yield_trkVeto += myLumi*scale1fb();
	  yGen_trkVeto++;

	  // Tau veto
	  if( !PassTauVeto() ) continue;
	  yield_tauVeto += myLumi*scale1fb();
	  yGen_tauVeto++;

	  // 4-jet requirement
	  if( ngoodjets() < 4 ) continue;
	  yield_4jets += myLumi*scale1fb();
	  yGen_4jets++;

	  // B-tag requirement
	  // if( ngoodbtags() < 1 ) continue;
	  int nbtags = 0;
	  for( uint idx=0; idx<ak4pfjets_CSV().size(); idx++ ) {
		if( ak4pfjets_pt().at(idx) <= 30. ) continue;
		if( fabs(ak4pfjets_eta().at(idx)) >= 2.4 ) continue;
		if( !ak4pfjets_loose_pfid().at(idx) ) continue;
		if(        TString(sampleName).Contains("stop")  && ak4pfjets_CSV().at(idx) < 0.814 ) continue; // Different CSV cut for phys14 signal samples
		else if( !(TString(sampleName).Contains("stop")) && ak4pfjets_CSV().at(idx) < 0.890 ) continue;
		nbtags++;
	  }
	  if( nbtags < 1 ) continue;
	  yield_1bjet += myLumi*scale1fb();
	  yGen_1bjet++;

	  // Baseline MET cut
	  if( pfmet() <= 200. ) continue;
	  yield_METcut += myLumi*scale1fb();
	  yGen_METcut++;

	  // MT cut
	  if( mt_met_lep() <= 150. ) continue;
	  yield_MTcut += myLumi*scale1fb();
	  yGen_MTcut++;

	  // Min delta-phi between MET and j1/j2
	  if( mindphi_met_j1_j2() <= 0.8 ) continue;
	  yield_dPhi += myLumi*scale1fb();
	  yGen_dPhi++;

	  // Chi^2 cut
	  // if( hadronic_top_chi2() >= 10. ) continue;
	  yield_chi2 += myLumi*scale1fb();
	  yGen_chi2++;


	  //////////////////////////////////////////////////////////
	  // Classify event based on number of leptons / neutrinos

	  int nPromptLeps = 0;
	  bool nusFromZ = false;

	  if( !TString(sampleName).Contains("stop") ) {

		for( uint idx=0; idx<genels_p4().size(); idx++ ) {
		  if( abs(genels_motherid().at(idx)) == 15 ) continue;  // Skip leptons from tau decay
		  if( !genels_fromHardProcessFinalState().at(idx) ) continue;
		  if( !genels_isLastCopy().at(idx) ) continue;
		  nPromptLeps++;
		}

		for( uint idx=0; idx<genmus_p4().size(); idx++ ) {
		  if( abs(genmus_motherid().at(idx)) == 15 ) continue;  // Skip leptons from tau decay
		  if( !genmus_fromHardProcessFinalState().at(idx) ) continue;
		  if( !genmus_isLastCopy().at(idx) ) continue;
		  nPromptLeps++;
		}

		for( uint idx=0; idx<gentaus_p4().size(); idx++ ) {
		  if( !gentaus_fromHardProcessDecayed().at(idx) ) continue;
		  if( !gentaus_isLastCopy().at(idx) ) continue;
		  nPromptLeps++;
		}

		vector<int> nuIDs, nuMoms;

		// Make lists of neutrino IDs and parent indices
		for( uint idx=0; idx<gennuels_p4().size(); idx++  ) {
		  if( !gennuels_fromHardProcessFinalState().at(idx) ) continue;
		  if( !gennuels_isLastCopy().at(idx)  ) continue;
		  nuIDs.push_back( gennuels_id().at(idx) );
		  nuMoms.push_back( gennuels_motheridx().at(idx) );
		}
		for( uint idx=0; idx<gennumus_p4().size(); idx++  ) {
		  if( !gennumus_fromHardProcessFinalState().at(idx) ) continue;
		  if( !gennumus_isLastCopy().at(idx)  ) continue;
		  nuIDs.push_back( gennumus_id().at(idx) );
		  nuMoms.push_back( gennumus_motheridx().at(idx) );
		}
		for( uint idx=0; idx<gennutaus_p4().size(); idx++ ) {
		  if( !gennutaus_fromHardProcessFinalState().at(idx) ) continue;
		  if( !gennutaus_isLastCopy().at(idx) ) continue;
		  nuIDs.push_back( gennutaus_id().at(idx) );
		  nuMoms.push_back( gennutaus_motheridx().at(idx) );
		}


		// Hunt for two neutrinos with opposite PDG IDs and the same mother idx
		uint upperBound = nuIDs.size()<2 ?  0  :  nuIDs.size()-1;
		for( uint x=0; x<upperBound; x++ ) {
		  for( uint y=x+1; y<nuIDs.size(); y++ ) {
			if( (nuIDs.at(x) + nuIDs.at(y) == 0) && (nuMoms.at(x) == nuMoms.at(y) ) ) nusFromZ = true;
		  }
		}

	  } // End "if sampleName doesn't contain 'stop'"


	  int category = -99;
	  if(      nusFromZ==true   ) category = 3;   // Z to nu nu
	  else if( nPromptLeps >= 2 ) category = 2;   // 2 or more leptons
	  else if( nPromptLeps == 1 ) category = 1;   // 1 lepton
	  else                        category = 4;   // Other

	  ///////////////////////////////////////////
	  // MET/MT2W cuts and histo filling

	  // If the event passes the SR cuts, store which background type this event is, and fill histograms
	  for( int i=0; i<9; i++ ) {
		if( pfmet() < met_min[i] || pfmet() >= met_max[i] ) continue;
		if( MT2W() < mt2w_min[i] || MT2W() >= mt2w_max[i] ) continue;

		h_bgtype[i]->Fill( category,                    myLumi*scale1fb() );

		h_mt[i]->Fill(      mt_met_lep(), 				myLumi*scale1fb() );
		h_met[i]->Fill(     pfmet(),					myLumi*scale1fb() );
		h_mt2w[i]->Fill(	MT2W(),   					myLumi*scale1fb() );
		h_chi2[i]->Fill(	hadronic_top_chi2(),		myLumi*scale1fb() );
		h_htratio[i]->Fill( ak4_htratiom(),				myLumi*scale1fb() );
		h_mindphi[i]->Fill( mindphi_met_j1_j2() ,		myLumi*scale1fb() );
		h_ptb1[i]->Fill(	ak4pfjets_leadMEDbjet_pt(),	myLumi*scale1fb() );
		h_drlb1[i]->Fill(   dR_lep_leadb(),				myLumi*scale1fb() );
		h_ptlep[i]->Fill(   lep1_pt(),					myLumi*scale1fb() );
		h_metht[i]->Fill(   MET_over_sqrtHT(),			myLumi*scale1fb() );
		h_dphilw[i]->Fill(  dphi_Wlep(),				myLumi*scale1fb() );
		h_njets[i]->Fill(   ngoodjets(),                myLumi*scale1fb() );
		h_nbtags[i]->Fill(  nbtags,                     myLumi*scale1fb() );

		h_sigRegion->Fill( float(i+1),                  myLumi*scale1fb() );
	  }

	  // ---------------------------------------------------------------------------------------------------//
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////
    } //End of loop over events
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }

  cout << "Cutflow yields:                        (yield)  (gen evts)" << endl;

  printf("Total number of events:             %10.2f %9i\n", yield_total	, yGen_total		);
  printf("Events with 1st vertex good:        %10.2f %9i\n", yield_vtx		, yGen_vtx			);
  printf("Events with at least 1 good lepton: %10.2f %9i\n", yield_1goodlep	, yGen_1goodlep		);
  printf("Events passing second lepton veto:  %10.2f %9i\n", yield_2lepveto	, yGen_2lepveto		);
  printf("Events passing lepton selection:    %10.2f %9i\n", yield_lepSel	, yGen_lepSel		);
  printf("Events passing track veto:          %10.2f %9i\n", yield_trkVeto	, yGen_trkVeto		);
  printf("Events passing tau veto:            %10.2f %9i\n", yield_tauVeto	, yGen_tauVeto		);
  printf("Events with at least 4 jets:        %10.2f %9i\n", yield_4jets	, yGen_4jets		);
  printf("Events with at least 1 b-tag:       %10.2f %9i\n", yield_1bjet	, yGen_1bjet		);
  printf("Events with MET > 200 GeV:          %10.2f %9i\n", yield_METcut	, yGen_METcut		);
  printf("Events with MT > 150 GeV:           %10.2f %9i\n", yield_MTcut	, yGen_MTcut		);
  printf("Events with min dPhi > 0.8:         %10.2f %9i\n", yield_dPhi		, yGen_dPhi			);
  // printf("Events with chi2 < 10:              %10.2f %9i\n", yield_chi2 	, yGen_chi2 		);
  printf("Yield after preselection:           %10.2f %9i\n", yield_chi2		, yGen_chi2			);

  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  // Zero negative values in each signal region
  for( int j=0; j<9; j++ ) {
	bool negsFound = false;

	// First zero any decay modes with negative yields
	for( int k=1; k<= h_bgtype[j]->GetNbinsX(); k++ ) {
	  if( h_bgtype[j]->GetBinContent(k) < 0.0 ) {
		h_bgtype[j]->SetBinContent(k, 0.);
		h_bgtype[j]->SetBinError(k, 0.);
		negsFound = true;
	  }
	}
	// If any negative yields were found in any decay mode, recalculate the total yield
	if( negsFound ) {
	  double newYield, newErr;
	  newYield = h_bgtype[j]->IntegralAndError( 0, -1, newErr );
	  h_sigRegion->SetBinContent(j+1, newYield);
	  h_sigRegion->SetBinError(j+1, newErr);
	}
  }

  // Store histograms and clean them up
  TFile* plotfile = new TFile("plots.root", "UPDATE");
  plotfile->cd();

  for( int j=0; j<9; j++ ) {
	h_bgtype[j]->Write();
	h_mt[j]->Write();
	h_met[j]->Write();
	h_mt2w[j]->Write();
	h_chi2[j]->Write();
	h_htratio[j]->Write();
	h_mindphi[j]->Write();
	h_ptb1[j]->Write();
	h_drlb1[j]->Write();
	h_ptlep[j]->Write();
	h_metht[j]->Write();
	h_dphilw[j]->Write();
	h_njets[j]->Write();
	h_nbtags[j]->Write();
  }
  h_sigRegion->Write();

  plotfile->Close();

  for( int j=0; j<9; j++ ) {
	h_bgtype[j]->Delete();
	h_mt[j]->Delete();
	h_met[j]->Delete();
	h_mt2w[j]->Delete();
	h_chi2[j]->Delete();
	h_htratio[j]->Delete();
	h_mindphi[j]->Delete();
	h_ptb1[j]->Delete();
	h_drlb1[j]->Delete();
	h_ptlep[j]->Delete();
	h_metht[j]->Delete();
	h_dphilw[j]->Delete();
	h_njets[j]->Delete();
	h_nbtags[j]->Delete();
  }
  h_sigRegion->Delete();

  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
}
