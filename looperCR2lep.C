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
#include "TH2F.h"
#include "TH3D.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TVector3.h"
#include "Math/VectorUtil.h"

// CMS3
#include "CMS3.h"

// Custom
#include "analysis.h"
#include "sample.h"
#include "sfManager.h"

using namespace std;
using namespace tas;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;


int looperCR2lep( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true) {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Setup
  TChain *chain = mySample->GetChain();
  TString sampleName = mySample->GetLabel();
  const double myLumi = myAnalysis->GetLumi();
  const int nSigRegs = myAnalysis->GetSigRegionsAll().size();
  bool isFastsim = mySample->IsSignal();
  cout << "\nSample: " << sampleName.Data() << endl;

  /////////////////////////////////////////////////////////
  // Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  TH1::SetDefaultSumw2();


  TH1D* h_bgtype[nSigRegs];

  TH1D *h_mt[nSigRegs];
  TH1D *h_met[nSigRegs];
  TH1D *h_mt2w[nSigRegs];
  TH1D *h_chi2[nSigRegs];
  TH1D *h_htratio[nSigRegs];
  TH1D *h_mindphi[nSigRegs];
  TH1D *h_ptb1[nSigRegs];
  TH1D *h_drlb1[nSigRegs];
  TH1D *h_ptlep[nSigRegs];
  TH1D *h_metht[nSigRegs];
  TH1D *h_dphilw[nSigRegs];
  TH1D *h_njets[nSigRegs];
  TH1D *h_nbtags[nSigRegs];

  double met_min[nSigRegs]   = {250., 350., 250., 350.,   250., 325., 250., 350., 450., 250.};
  double met_max[nSigRegs]   = {350., 99999., 350., 99999., 325., 99999., 350., 450., 99999., 99999.};
  double mt2w_min[nSigRegs]  = {  0.,   0., 200., 200.,     0.,   0., 200., 200., 200., 0.};
  double mt2w_max[nSigRegs]  = {99999., 99999., 99999., 99999., 200., 200., 99999., 99999., 99999., 99999,};
  double njets_min[nSigRegs] = {   2,    2,    3,    3,      4,    4,    4,    4,    4, 2};
  double njets_max[nSigRegs] = {   2,    2,    3,    3,    999,  999,  999,  999,  999, 999};
  string regNames[nSigRegs] = {"compr250CR", "compr350CR", "boost250CR", "boost350CR", "low250CR", "low325CR", "high250CR", "high350CR", "high450CR", "inclusiveCR"};


  for( int i=0; i<nSigRegs; i++ ) {

	h_bgtype[i]   = new TH1D( Form( "bkgtype_%s_%s" , sampleName.Data(), regNames[i].c_str()), "Yield by background type",  4, 0.5, 4.5);
	h_mt[i]       = new TH1D( Form( "mt_%s_%s"      , sampleName.Data(), regNames[i].c_str()), "Transverse mass",			80, 0, 800);
	h_met[i]      = new TH1D( Form( "met_%s_%s"     , sampleName.Data(), regNames[i].c_str()), "MET",						40, 0, 1000);
	h_mt2w[i]     = new TH1D( Form( "mt2w_%s_%s"    , sampleName.Data(), regNames[i].c_str()), "MT2W",						50, 0, 500);
	h_chi2[i]     = new TH1D( Form( "chi2_%s_%s"    , sampleName.Data(), regNames[i].c_str()), "Hadronic #chi^{2}", 		50, 0, 15);
	h_htratio[i]  = new TH1D( Form( "htratio_%s_%s" , sampleName.Data(), regNames[i].c_str()), "H_{T} ratio",				50, 0, 1);
	h_mindphi[i]  = new TH1D( Form( "mindphi_%s_%s" , sampleName.Data(), regNames[i].c_str()), "min #Delta#phi(j12,MET)",	50, 0, 4);
	h_ptb1[i]     = new TH1D( Form( "ptb1_%s_%s"    , sampleName.Data(), regNames[i].c_str()), "p_{T} (b1)",				100, 0, 500);
	h_drlb1[i]    = new TH1D( Form( "drlb1_%s_%s"   , sampleName.Data(), regNames[i].c_str()), "#DeltaR (lep, b1)", 		50, 0, 5);
	h_ptlep[i]    = new TH1D( Form( "ptlep_%s_%s"   , sampleName.Data(), regNames[i].c_str()), "p_{T} (lep)",				100, 0, 500);
	h_metht[i]    = new TH1D( Form( "metht_%s_%s"   , sampleName.Data(), regNames[i].c_str()), "MET/sqrt(HT)",				50, 0, 100);
	h_dphilw[i]   = new TH1D( Form( "dphilw_%s_%s"  , sampleName.Data(), regNames[i].c_str()), "#Delta#phi (lep,W)",		50, 0, 3.5);
	h_njets[i]    = new TH1D( Form( "njets_%s_%s"   , sampleName.Data(), regNames[i].c_str()), "Number of jets",            16, -0.5, 15.5);
	h_nbtags[i]   = new TH1D( Form( "nbtags_%s_%s"  , sampleName.Data(), regNames[i].c_str()), "Number of b-tags",          7, -0.5, 6.5);

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

  TH1D *h_yields = new TH1D( Form("srYields_%s", sampleName.Data()), "Yield by signal region", nSigRegs, 0.5, float(nSigRegs)+0.5);
  for( int i=0; i<nSigRegs; i++ ) h_yields->GetXaxis()->SetBinLabel( i+1, regNames[i].c_str() );
  h_yields->SetDirectory(rootdir);


  float yield_total = 0;
  float yield_vtx = 0;
  float yield_1goodlep = 0;
  // float yield_lepSel = 0;
  float yield_2lepveto = 0;
  float yield_trkVeto = 0;
  float yield_2lepCR = 0;
  float yield_tauVeto = 0;
  float yield_njets = 0;
  float yield_1bjet = 0;
  float yield_METcut = 0;
  float yield_MTcut = 0;
  float yield_dPhi = 0;
  float yield_chi2 = 0;

  int yGen_total = 0;
  int yGen_vtx = 0;
  int yGen_1goodlep = 0;
  // int yGen_lepSel = 0;
  int yGen_2lepveto = 0;
  int yGen_trkVeto = 0;
  int yGen_tauVeto = 0;
  int yGen_2lepCR = 0;
  int yGen_njets = 0;
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
    TFile file( currentFile->GetTitle() );
    TTree *tree = (TTree*)file.Get("t");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    cms3.Init(tree);

	// Initialize scale factor manager
	sfManager mySFs( isFastsim, ".", (TH1D*)file.Get( "h_counter" ) );
    
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
	  // if(      sampleName == "tt2l"  && genlepsfromtop() != 2 ) continue;  //Require 2 leps from top in "tt2l" events
	  // else if( sampleName == "tt1l"  && genlepsfromtop() != 1 ) continue;  //Require 1 lep from top in "tt1l" events


	  /////////////////////////////////
	  // Set event weight

	  double evtWeight = 1.;

	  if( is_data() || mySample->IsData() ) evtWeight = 1.;
	  else if( mySample->IsSignal() ) {
		TH2F* hNEvts = (TH2F*)file.Get("histNEvts");
		TH3D* hCounter = (TH3D*)file.Get("h_counterSMS");
		double nEvtsSample = hNEvts->GetBinContent( hNEvts->FindBin( mass_stop(), mass_lsp() ) );
		int binx = hCounter->GetXaxis()->FindBin( mass_stop() );
		int biny = hCounter->GetYaxis()->FindBin( mass_lsp()  );
		double bTagSumWeights = hCounter->GetBinContent( binx, biny, 14 );
		mySFs.SetBtagNorm( nEvtsSample / bTagSumWeights );
		evtWeight = myLumi * 1000. * xsec() / nEvtsSample;
	  }
	  else evtWeight = myLumi * scale1fb();


	  // Count the number of events processed
	  yield_total += evtWeight;
	  yGen_total++;

	  // First vertex must be good
	  if( firstGoodVtxIdx() != 0 ) continue;
	  yield_vtx += evtWeight;
	  yGen_vtx++;

	  // Must have at least 1 good lepton
	  if( ngoodleps() < 1 ) continue;
	  if(      !is_data() && lep1_is_el() ) evtWeight *= mySFs.GetSF_el( lep1_pt(), lep1_eta() );
	  else if( !is_data() && lep1_is_mu() ) evtWeight *= mySFs.GetSF_mu( lep1_pt(), lep1_eta() );
	  yield_1goodlep += evtWeight;
	  yGen_1goodlep++;

	  // Lep 1 must pass lepton selections
	  // These aren't necessary unless syncing with John
	  // if( lep1_is_el() ) {
	  // 	if( lep1_pt() < 20. ) continue;
	  // 	if( fabs(lep1_eta()) > 1.442 ) continue;  // It's 1.4442 in the babymaker, but John uses 1.442
	  // 	if( !lep1_is_eleid_medium() ) continue;
	  // }
	  // else if( lep1_is_mu() ) {
	  // 	if( lep1_pt() < 20. ) continue;
	  // 	if( fabs(lep1_eta()) > 2.4 ) continue;
	  // 	if( !lep1_is_muoid_tight() ) continue;
	  // }
	  // yield_lepSel += evtWeight;
	  // yGen_lepSel++;


	  ///////////////////
	  // Make 2-lepton CR

	  int countGoodLeps = 0;

	  // Second lepton veto
	  // if( nvetoleps() > 1 && ROOT::Math::VectorUtil::DeltaR( lep1_p4(), lep2_p4() ) > 0.01 ) continue;
	  // yield_2lepveto += evtWeight;
	  // yGen_2lepveto++;

	  // Count the number of veto leptons, but subtract one if lep1 and lep2 overlap
	  countGoodLeps += nvetoleps();

	  if( nvetoleps() == 2 && ROOT::Math::VectorUtil::DeltaR( lep1_p4(), lep2_p4() ) < 0.01 ) countGoodLeps--;
	  else if( nvetoleps() >= 2 && lep2_pt() < 10. ) countGoodLeps = 1;

	  if( !is_data() && nvetoleps() >= 2 && lep2_pt() >= 10. ) {
		if(      lep2_is_el() ) evtWeight *= mySFs.GetSF_elVeto( lep2_pt(), lep2_eta() );
		else if( lep2_is_mu() ) evtWeight *= mySFs.GetSF_muVeto( lep2_pt(), lep2_eta() );
	  }

	  if( countGoodLeps > 1 ) {
		  yield_2lepveto += evtWeight;
		  yGen_2lepveto++;
	  }

	  // Track veto
	  // if( !PassTrackVeto_v3() ) continue;
	  // yield_trkVeto += evtWeight;
	  // yGen_trkVeto++;

	  // If we fail the track veto, count another good lepton
	  if( !PassTrackVeto_v3() ) {
		countGoodLeps++;
		yield_trkVeto += evtWeight;
		yGen_trkVeto++;
	  }

	  // Tau veto
	  // if( !PassTauVeto() ) continue;
	  // yield_tauVeto += evtWeight;
	  // yGen_tauVeto++;

	  if( !PassTauVeto() ) {
		countGoodLeps++;
		yield_tauVeto += evtWeight;
		yGen_tauVeto++;
	  }

	  if( countGoodLeps < 2 ) continue;
	  yield_2lepCR += evtWeight;
	  yGen_2lepCR++;



	  ////////////////////
	  ////////////////////

	  // N-jet requirement
	  if( ngoodjets() < 2 ) continue;
	  if( !is_data() ) {
		for( uint i=0; i<ak4pfjets_p4().size(); i++ ) {
		  evtWeight *= mySFs.GetSF_btag( ak4pfjets_pt().at(i),
										 ak4pfjets_eta().at(i),
										 ak4pfjets_hadron_flavor().at(i),
										 ak4pfjets_CSV().at(i) );
		}
	  }
	  yield_njets += evtWeight;
	  yGen_njets++;

	  // B-tag requirement
	  if( ngoodbtags() < 1 ) continue;
	  yield_1bjet += evtWeight;
	  yGen_1bjet++;

	  // Baseline MET cut
	  if( pfmet() <= 250. ) continue;
	  yield_METcut += evtWeight;
	  yGen_METcut++;

	  // MT cut
	  if( mt_met_lep() <= 150. ) continue;
	  yield_MTcut += evtWeight;
	  yGen_MTcut++;

	  // Min delta-phi between MET and j1/j2
	  if( mindphi_met_j1_j2() <= 0.8 ) continue;
	  yield_dPhi += evtWeight;
	  yGen_dPhi++;

	  // Chi^2 cut
	  // if( hadronic_top_chi2() >= 10. ) continue;
	  yield_chi2 += evtWeight;
	  yGen_chi2++;


	  //////////////////////////////////////////////////////////
	  // Classify event based on number of leptons / neutrinos

	  int category = -99;
	  if(   isZtoNuNu() ) category = 3;   // Z to nu nu
	  else if( is2lep() ) category = 2;   // 2 or more leptons
	  else if( is1lep() ) category = 1;   // 1 lepton
	  else                category = 4;   // Other

	  ///////////////////////////////////////////
	  // MET/MT2W cuts and histo filling

	  // If the event passes the SR cuts, store which background type this event is, and fill histograms
	  for( int i=0; i<nSigRegs; i++ ) {
		if( pfmet() < met_min[i] || pfmet() >= met_max[i] ) continue;
		if( MT2W() < mt2w_min[i] || MT2W() >= mt2w_max[i] ) continue;
		if( ngoodjets() < njets_min[i] || ngoodjets() > njets_max[i] ) continue;
		if( (regNames[i] == "compr250" || regNames[i] == "compr350") && topnessMod() <= 6.4 ) continue;

		h_bgtype[i]->Fill( category,                    evtWeight );

		h_mt[i]->Fill(      mt_met_lep(), 				evtWeight );
		h_met[i]->Fill(     pfmet(),					evtWeight );
		h_mt2w[i]->Fill(	MT2W(),   					evtWeight );
		h_chi2[i]->Fill(	hadronic_top_chi2(),		evtWeight );
		h_htratio[i]->Fill( ak4_htratiom(),				evtWeight );
		h_mindphi[i]->Fill( mindphi_met_j1_j2() ,		evtWeight );
		h_ptb1[i]->Fill(	ak4pfjets_leadMEDbjet_pt(),	evtWeight );
		h_drlb1[i]->Fill(   dR_lep_leadb(),				evtWeight );
		h_ptlep[i]->Fill(   lep1_pt(),					evtWeight );
		h_metht[i]->Fill(   MET_over_sqrtHT(),			evtWeight );
		h_dphilw[i]->Fill(  dphi_Wlep(),				evtWeight );
		h_njets[i]->Fill(   ngoodjets(),                evtWeight );
		h_nbtags[i]->Fill(  ngoodbtags(),               evtWeight );

		h_yields->Fill(     float(i+1),                 evtWeight );
	  }

	  // ---------------------------------------------------------------------------------------------------//
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////
    } //End of loop over events
  
    // Clean Up
    delete tree;
    file.Close();
  }

  cout << "Cutflow yields:                        (yield)  (gen evts)" << endl;

  printf("Total number of events:             %10.2f %9i\n", yield_total	, yGen_total		);
  printf("Events with 1st vertex good:        %10.2f %9i\n", yield_vtx		, yGen_vtx			);
  printf("Events with at least 1 good lepton: %10.2f %9i\n", yield_1goodlep	, yGen_1goodlep		);
  // printf("Events passing lepton selection:    %10.2f %9i\n", yield_lepSel	, yGen_lepSel		);

  printf("\nEvents passing 2-lep requirement:   %10.2f %9i\n", yield_2lepCR   , yGen_2lepCR       );
  printf("   Events with veto lepton:         %10.2f %9i\n", yield_2lepveto	, yGen_2lepveto		);
  printf("   Events with isolated track:      %10.2f %9i\n", yield_trkVeto	, yGen_trkVeto		);
  printf("   Events with identified tau:      %10.2f %9i\n\n", yield_tauVeto	, yGen_tauVeto		);

  printf("Events with at least 2 jets:        %10.2f %9i\n", yield_njets	, yGen_njets		);
  printf("Events with at least 1 b-tag:       %10.2f %9i\n", yield_1bjet	, yGen_1bjet		);
  printf("Events with MET > 250 GeV:          %10.2f %9i\n", yield_METcut	, yGen_METcut		);
  printf("Events with MT > 150 GeV:           %10.2f %9i\n", yield_MTcut	, yGen_MTcut		);
  printf("Events with min dPhi > 0.8:         %10.2f %9i\n", yield_dPhi		, yGen_dPhi			);
  // printf("Events with chi2 < 10:              %10.2f %9i\n", yield_chi2 	, yGen_chi2 		);
  printf("Yield after preselection:           %10.2f %9i\n", yield_chi2		, yGen_chi2			);

  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  // Zero negative values in each signal region
  for( int j=0; j<nSigRegs; j++ ) {
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
	  h_yields->SetBinContent(j+1, newYield);
	  h_yields->SetBinError(j+1, newErr);
	}
  }

  // Store histograms and clean them up
  TFile* plotfile = new TFile( myAnalysis->GetFileName(), "UPDATE");
  plotfile->cd();

  for( int j=0; j<nSigRegs; j++ ) {
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
  h_yields->Write();

  plotfile->Close();

  for( int j=0; j<nSigRegs; j++ ) {
	delete h_bgtype[j];
	delete h_mt[j];
	delete h_met[j];
	delete h_mt2w[j];
	delete h_chi2[j];
	delete h_htratio[j];
	delete h_mindphi[j];
	delete h_ptb1[j];
	delete h_drlb1[j];
	delete h_ptlep[j];
	delete h_metht[j];
	delete h_dphilw[j];
	delete h_njets[j];
	delete h_nbtags[j];
  }
  delete h_yields;

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
