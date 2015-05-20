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

// CMS3
#include "CMS3.cc"

using namespace std;
using namespace tas;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;


int ScanChain( TChain* chain, int nLepRequired = -1, bool fast = true, int nEvents = -1 /*, string skimFilePrefix = "test"*/) {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  /////////////////////////////////////////////////////////////////////
  // Example Histograms
  // TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  // TH1F *samplehisto = new TH1F("samplehisto", "Example histogram", 200,0,200);
  // samplehisto->SetDirectory(rootdir);

  int nEvt_total = 0;
  int nEvt_vtx = 0;
  int nEvt_goodlep = 0;
  int nEvt_1goodlep = 0;
  int nEvt_lepSel = 0;
  int nEvt_lepIso = 0;
  int nEvt_trkVeto = 0;
  int nEvt_tauVeto = 0;
  int nEvt_4jets = 0;
  int nEvt_1bjet = 0;
  int nEvt_METcut = 0;
  int nEvt_MTcut = 0;
  int nEvt_dPhi = 0;
  int nEvt_chi2 = 0;

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
	  if( nLepRequired > 0 && nGenLep != nLepRequired ) continue;


	  // Count the number of events
	  nEvt_total++;

	  // Must have >=1 good vertex
	  if( nvtxs() < 1 ) continue;
	  nEvt_vtx++;

	  // Must have at least 1 good lepton
	  if( ngoodlep() < 1 ) continue;
	  nEvt_goodlep++;

	  // Must have EXACTLY 1 good lepton;
	  if( ngoodlep() > 1 ) continue;
	  nEvt_1goodlep++;

	  // Must pass lepton selections
	  if( lep1_is_el() ) {
		if( lep1_pt() < 30. ) continue;
		if( fabs(lep1_eta()) >= 1.442 ) continue;
		if( !lep1_is_eleid_medium() ) continue;
	  }
	  else if( lep1_is_mu() ) {
		if( lep1_pt() < 25. ) continue;
		if( fabs(lep1_eta()) >= 2.1 ) continue;
		if( !lep1_is_muoid_tight() ) continue;
	  }
	  nEvt_lepSel++;

	  // Lepton isolation requirement
	  // pt_sum < min( 5GeV, 0.15*lep_pT )
	  double lepiso = lep1_relIso03DB() * lep1_pt();
	  double mymin = fmin( 5, 0.15*lep1_pt()  );
	  if( lepiso > mymin ) continue;
	  nEvt_lepIso++;

	  // Track veto
	  if( !PassTrackVeto() ) continue;
	  nEvt_trkVeto++;

	  // Tau veto
	  bool foundtau = false;
	  for( uint i=0; i<tau_isVetoTau().size(); i++ ) {
		if( tau_isVetoTau().at(i) ) foundtau=true;
	  }
	  if( foundtau ) continue;
	  nEvt_tauVeto++;

	  // 4-jet requirement
	  int ngoodjets = 0;
	  int ngoodbjets = 0;
	  for( uint i=0; i<ak4pfjets_loose_pfid().size(); i++ ) {
		if( ak4pfjets_p4().at(i).pt() <= 30. ) continue;
		if( fabs(ak4pfjets_p4().at(i).eta()) >= 2.4 ) continue;
		if( !ak4pfjets_loose_pfid().at(i) ) continue;
		// if( !ak4pfjets_loose_puid().at(i) ) continue;
		ngoodjets++;
		if( !ak4pfjets_passMEDbtag().at(i) ) continue;
		ngoodbjets++;
	  }
	  if( ngoodjets < 4 ) continue;
	  nEvt_4jets++;

	  // B-tag requirement
	  if( ngoodbjets < 1 ) continue;
	  nEvt_1bjet++;

	  // MET requirement
	  if( pfmet() <= 100. ) continue;
	  nEvt_METcut++;

	  // MT requirement
	  if( MT_MET_lep1() <= 120. ) continue;
	  nEvt_MTcut++;

	  // Delta-Phi cut
	  if( mindphi_met_j1_j2() <= 0.8 ) continue;
	  nEvt_dPhi++;

	  // Chi^2 cut
	  if( chi2() > 5. ) continue;
	  nEvt_chi2++;





	  // ---------------------------------------------------------------------------------------------------//
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }

  cout << "Cutflow yields:" << endl;

  printf("Total number of events:             %8d\n", nEvt_total );
  printf("Events with at least 1 good vertex: %8d\n", nEvt_vtx );
  printf("Events with at least 1 good lepton: %8d\n", nEvt_goodlep );
  printf("Events with exactly 1 good lepton:  %8d\n", nEvt_1goodlep );
  printf("Events passing lepton selection:    %8d\n", nEvt_lepSel );
  printf("Events passing lepton isolation:    %8d\n", nEvt_lepIso );
  printf("Events passing track veto:          %8d\n", nEvt_trkVeto );
  printf("Events passing tau veto:            %8d\n", nEvt_tauVeto );
  printf("Events with at least 4 jets:        %8d\n", nEvt_4jets );
  printf("Events with at least 1 b-tag:       %8d\n", nEvt_1bjet );
  printf("Events with MET > 100 GeV:          %8d\n", nEvt_METcut );
  printf("Events with MT > 120 GeV:           %8d\n", nEvt_MTcut );
  printf("Events with min dPhi > 0.8:         %8d\n", nEvt_dPhi );
  printf("Events with chi2 < 5:               %8d\n", nEvt_chi2 );

  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  
  // Example Histograms
  // samplehisto->Draw();
  
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
