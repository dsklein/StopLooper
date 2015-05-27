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

// CMS3
#include "CMS3.cc"

using namespace std;
using namespace tas;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;


int ScanChain( TChain* chain, string samplename = "noname", int nLepRequired = -1, bool fast = true, int nEvents = -1 /*, string skimFilePrefix = "test"*/) {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  /////////////////////////////////////////////////////////////////////
  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  TH1::SetDefaultSumw2();

  TH1D *h_mt[8];
  TH1D *h_met[8];
  TH1D *h_mt2w[8];
  TH1D *h_chi2[8];
  TH1D *h_htratio[8];
  TH1D *h_mindphi[8];
  TH1D *h_ptb1[8];
  TH1D *h_drlb1[8];
  TH1D *h_ptlep[8];
  TH1D *h_metht[8];
  TH1D *h_dphilw[8];

  int i = 0;

  for( float cut_mt2w : {0., 200.} ) {
	for( float cut_met : {150., 200., 250., 300.} ) {

	  TString s_mt      = "mt_"      + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_met     = "met_"     + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_mt2w    = "mt2w_"    + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_chi2    = "chi2_"    + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_htratio = "htratio_" + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_mindphi = "mindphi_" + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_ptb1    = "ptb1_"    + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_drlb1   = "drlb1_"   + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_ptlep   = "ptlep_"   + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_metht   = "metht_"   + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);
	  TString s_dphilw  = "dphilw_"  + samplename + Form("_%.f_%3.f", cut_mt2w, cut_met);

	  h_mt[i]       = new TH1D( s_mt.Data(),		"Transverse mass",			50, 100, 600);
	  h_met[i]      = new TH1D( s_met.Data(),		"MET",						50, 0, 1000);
	  h_mt2w[i]     = new TH1D( s_mt2w.Data(),		"MT2W",						50, 0, 500);
	  h_chi2[i]     = new TH1D( s_chi2.Data(),		"Hadronic #chi^{2}",		50, 0, 6);
	  h_htratio[i]  = new TH1D( s_htratio.Data(),	"H_{T} ratio",				50, 0, 1);
	  h_mindphi[i]  = new TH1D( s_mindphi.Data(),	"min #Delta#phi(j12,MET)",	50, 0, 4);
	  h_ptb1[i]     = new TH1D( s_ptb1.Data(),		"p_{T} (b1)",				50, 0, 800);
	  h_drlb1[i]    = new TH1D( s_drlb1.Data(),		"#DeltaR (lep, b1)",		50, 0, 5);
	  h_ptlep[i]    = new TH1D( s_ptlep.Data(),		"p_{T} (lep)",				50, 0, 800);
	  h_metht[i]    = new TH1D( s_metht.Data(),		"MET/sqrt(HT)",				50, 0, 100);
	  h_dphilw[i]   = new TH1D( s_dphilw.Data(),	"#Delta#phi (lep,W)",		50, 0, 3.5);

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

	  i++;
	}
  }

  TString yieldName = "yield_" + samplename;
  TH1D *h_yield = new TH1D( yieldName.Data(), "Yield", 9, -0.5, 8.5);
  h_yield->SetDirectory(rootdir);

  float yield_total = 0;
  float yield_vtx = 0;
  float yield_goodlep = 0;
  float yield_1goodlep = 0;
  float yield_lepSel = 0;
  float yield_lepIso = 0;
  float yield_trkVeto = 0;
  float yield_tauVeto = 0;
  float yield_4jets = 0;
  float yield_1bjet = 0;
  float yield_METcut = 0;
  float yield_MTcut = 0;
  float yield_dPhi = 0;
  float yield_chi2 = 0;


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


	  // First things first:
	  // If this is a top event, let's make sure it has the right number of gen leptons in it!
	  int nGenLep = gen_nfromtels_() + gen_nfromtmus_() + gen_nfromttaus_();
	  if( nLepRequired >= 0 && nGenLep != nLepRequired ) continue;


	  // Count the number of events
	  yield_total += 10.*scale1fb();

	  // Must have >=1 good vertex
	  if( nvtxs() < 1 ) continue;
	  yield_vtx += 10.*scale1fb();

	  // Must have at least 1 good lepton
	  if( ngoodlep() < 1 ) continue;
	  yield_goodlep += 10.*scale1fb();

	  // Must have EXACTLY 1 good lepton;
	  if( ngoodlep() > 1 ) continue;
	  yield_1goodlep += 10.*scale1fb();

	  // Must pass lepton selections
	  if( lep1_is_el() ) {
		if( lep1_pt() <= 40. ) continue;
		if( fabs(lep1_eta()) >= 2.1 ) continue;
		if( !lep1_is_phys14_medium_noIso() ) continue;
	  }
	  else if( lep1_is_mu() ) {
		if( lep1_pt() < 30. ) continue;
		if( fabs(lep1_eta()) >= 2.1 ) continue;
		// if( !lep1_is_muoid_tight() ) continue;   // This is meant to be medium, but there's no branch for medium in the babies right now
		if( fabs(lep1_d0()) >= 0.02 ) continue;  // Indara says this is buggy
		if( fabs(lep1_dz()) >= 0.1  ) continue;  // Indara says this is buggy
	  }
	  yield_lepSel += 10.*scale1fb();

	  // Lepton isolation requirement
	  if( lep1_miniRelIsoDB() >= 0.1 ) continue;
	  yield_lepIso += 10.*scale1fb();

	  // Track veto
	  if( !PassTrackVeto() ) continue;
	  yield_trkVeto += 10.*scale1fb();

	  // Tau veto
	  // bool foundtau = false;
	  // for( uint i=0; i<tau_isVetoTau().size(); i++ ) {
	  // 	if( tau_isVetoTau().at(i) ) foundtau=true;
	  // }
	  // if( foundtau ) continue;
	  if( !PassTauVeto() ) continue;
	  yield_tauVeto += 10.*scale1fb();

	  // 4-jet requirement
	  int ngoodjets = 0;
	  int ngoodbjets = 0;
	  for( uint j=0; j<ak4pfjets_loose_pfid().size(); j++ ) {
		if( ak4pfjets_p4().at(j).pt() <= 30. ) continue;
		if( fabs(ak4pfjets_p4().at(j).eta()) >= 2.4 ) continue;
		if( !ak4pfjets_loose_pfid().at(j) ) continue;
		// if( !ak4pfjets_loose_puid().at(j) ) continue;
		ngoodjets++;
		if( !ak4pfjets_passMEDbtag().at(j) ) continue;
		ngoodbjets++;
	  }
	  if( ngoodjets < 4 ) continue;
	  yield_4jets += 10.*scale1fb();

	  // B-tag requirement
	  if( ngoodbjets < 1 ) continue;
	  yield_1bjet += 10.*scale1fb();

	  // MET requirement
	  if( pfmet() <= 100. ) continue;
	  yield_METcut += 10.*scale1fb();

	  // MT requirement
	  if( MT_MET_lep1() <= 120. ) continue;
	  yield_MTcut += 10.*scale1fb();

	  // Delta-Phi cut
	  if( mindphi_met_j1_j2() <= 0.8 ) continue;
	  yield_dPhi += 10.*scale1fb();

	  // Chi^2 cut
	  if( chi2() >= 5. ) continue;
	  yield_chi2 += 10.*scale1fb();
	  h_yield->Fill( 0., 10.*scale1fb() );

	  //////////////////////////////////////////////////////
	  // Calculate a few variables that aren't in the babies

	  // MET / sqrt(HT)
	  double metSqHT = pfmet() / sqrt( ak4_HT() );

	  // Delta phi( lep, W )  -- Indara line 1137
	  double met = pfmet();
	  double phi = pfmet_phi();
	  const TVector3 lepVec( lep1_p4().x(), lep1_p4().y(), lep1_p4().z() );
	  const TVector3 metVec( met*cos(phi), met*sin(phi), 0 );
	  const TVector3 wVec = lepVec + metVec;
	  double dPhiLepMet = fabs( wVec.DeltaPhi(lepVec) );


	  //////////////////////////////////////////////////////
	  // Signal region cuts, and histo filling

	  i = -1;

	  for( float cut_mt2w : {0., 200.} ) {
		for( float cut_met : {150., 200., 250., 300.} ) {

		  i++;
		  if( MT2W_lep1() < cut_mt2w ) continue;
		  if( pfmet() < cut_met ) continue;

		  // Fill histograms!
		  h_mt[i]->Fill(      MT_MET_lep1(),				10.*scale1fb() );
		  h_met[i]->Fill(     pfmet(),						10.*scale1fb() );
		  h_mt2w[i]->Fill(	  MT2W_lep1(),					10.*scale1fb() );
		  h_chi2[i]->Fill(	  chi2(),						10.*scale1fb() );
		  h_htratio[i]->Fill( ak4_htratiom(),				10.*scale1fb() );
		  h_mindphi[i]->Fill( mindphi_met_j1_j2() ,			10.*scale1fb() );
		  h_ptb1[i]->Fill(	  ak4pfjets_leadMEDbjet_pt(),	10.*scale1fb() );
		  h_drlb1[i]->Fill(   dR_lep1_leadb(),				10.*scale1fb() );
		  h_ptlep[i]->Fill(   lep1_pt(),					10.*scale1fb() );
		  h_metht[i]->Fill(   metSqHT,						10.*scale1fb() );  // Calculate this yourself
		  h_dphilw[i]->Fill(  dPhiLepMet,					10.*scale1fb() );  // Calculate this yourself

		  h_yield->Fill( float(i+1), 10.*scale1fb() );
		}
	  }


	  // ---------------------------------------------------------------------------------------------------//
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }

  // cout << "Cutflow yields:" << endl;

  // printf("Total number of events:             %10.2f\n", yield_total );
  // printf("Events with at least 1 good vertex: %10.2f\n", yield_vtx );
  // printf("Events with at least 1 good lepton: %10.2f\n", yield_goodlep );
  // printf("Events with exactly 1 good lepton:  %10.2f\n", yield_1goodlep );
  // printf("Events passing lepton selection:    %10.2f\n", yield_lepSel );
  // printf("Events passing lepton isolation:    %10.2f\n", yield_lepIso );
  // printf("Events passing track veto:          %10.2f\n", yield_trkVeto );
  // printf("Events passing tau veto:            %10.2f\n", yield_tauVeto );
  // printf("Events with at least 4 jets:        %10.2f\n", yield_4jets );
  // printf("Events with at least 1 b-tag:       %10.2f\n", yield_1bjet );
  // printf("Events with MET > 100 GeV:          %10.2f\n", yield_METcut );
  // printf("Events with MT > 120 GeV:           %10.2f\n", yield_MTcut );
  // printf("Events with min dPhi > 0.8:         %10.2f\n", yield_dPhi );
  // printf("Events with chi2 < 5:               %10.2f\n", yield_chi2 );
  printf("Yield after preselection:           %10.2f\n", yield_chi2 );

  i = 2;
  for( float cut_mt2w : {0., 200.} ) {
	for( float cut_met : {150., 200., 250., 300.} ) {

	  double yield  = h_yield->GetBinContent(i);
	  double yError = h_yield->GetBinError(i);
	  printf("MT2W cut: %3.f, MET cut: %3.f,  yield: %8.1f +/- %6.1f\n", cut_mt2w, cut_met, yield, yError);
	  i++;
	}
  }

  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  
  // Write Histograms
  TFile* outfile = new TFile("plots.root", "UPDATE");
  outfile->cd();

  for( int j=0; j<8; j++ ) {
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
  }
  h_yield->Write();

  outfile->Close();

  for( int j=0; j<8; j++ ) {
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
  }
  h_yield->Delete();

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
