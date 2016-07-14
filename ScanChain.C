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
#include "TH2D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TVector3.h"
#include "Math/VectorUtil.h"

// CMS3 and CORE
#include "CMS3.h"
#include "../StopAnalysis_74x/CORE/Tools/badEventFilter.cc"
#include "../StopAnalysis_74x/CORE/Tools/dorky/dorky.cc"

// Custom
#include "analysis.h"
#include "sample.h"
#include "sfManager.h"

using namespace std;
using namespace tas;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

// Global variables used in defining signal regions
extern bool j1_isBtag;
extern double j1pt;


int ScanChain( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true) {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Setup
  TChain *chain = mySample->GetChain();
  TString sampleName = mySample->GetLabel();
  const int nSigRegs = myAnalysis->GetSigRegionsAll().size();
  bool isFastsim = mySample->IsSignal();
  cout << "\nSample: " << sampleName.Data() << endl;

  /////////////////////////////////////////////////////////
  // Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  TH1::SetDefaultSumw2();


  TH1D* h_bgtype[nSigRegs];
  TH1D* h_evttype[nSigRegs];
  TH2D* h_sigyields[nSigRegs];

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
  TH1D *h_ptj1[nSigRegs];
  TH1D *h_j1btag[nSigRegs];


  vector<TString> regNames = myAnalysis->GetSigRegionLabelsAll();
  vector<sigRegion> sigRegions = myAnalysis->GetSigRegionsAll();

  for( int i=0; i<nSigRegs; i++ ) {

	h_bgtype[i]   = new TH1D( Form( "bkgtype_%s_%s" , sampleName.Data(), regNames.at(i).Data()), "Yield by background type",  5, 0.5, 5.5);
	h_evttype[i]= new TH1D( Form( "evttype_%s"      , regNames.at(i).Data()),                    "Yield by event type",       6, 0.5, 6.5);
	h_sigyields[i] = new TH2D( Form( "sigyields_%s", regNames.at(i).Data()), "Signal yields by mass point", 37, 87.5, 1012.5, 21, -12.5, 512.5 );
	h_mt[i]       = new TH1D( Form( "mt_%s_%s"      , sampleName.Data(), regNames.at(i).Data()), "Transverse mass",			80, 0, 800);
	h_met[i]      = new TH1D( Form( "met_%s_%s"     , sampleName.Data(), regNames.at(i).Data()), "MET",						40, 0, 1000);
	h_mt2w[i]     = new TH1D( Form( "mt2w_%s_%s"    , sampleName.Data(), regNames.at(i).Data()), "MT2W",						50, 0, 500);
	h_chi2[i]     = new TH1D( Form( "chi2_%s_%s"    , sampleName.Data(), regNames.at(i).Data()), "Hadronic #chi^{2}", 		50, 0, 15);
	h_htratio[i]  = new TH1D( Form( "htratio_%s_%s" , sampleName.Data(), regNames.at(i).Data()), "H_{T} ratio",				50, 0, 1);
	h_mindphi[i]  = new TH1D( Form( "mindphi_%s_%s" , sampleName.Data(), regNames.at(i).Data()), "min #Delta#phi(j12,MET)",	50, 0, 4);
	h_ptb1[i]     = new TH1D( Form( "ptb1_%s_%s"    , sampleName.Data(), regNames.at(i).Data()), "p_{T} (b1)",				100, 0, 500);
	h_drlb1[i]    = new TH1D( Form( "drlb1_%s_%s"   , sampleName.Data(), regNames.at(i).Data()), "#DeltaR (lep, b1)", 		50, 0, 5);
	h_ptlep[i]    = new TH1D( Form( "ptlep_%s_%s"   , sampleName.Data(), regNames.at(i).Data()), "p_{T} (lep)",				100, 0, 500);
	h_metht[i]    = new TH1D( Form( "metht_%s_%s"   , sampleName.Data(), regNames.at(i).Data()), "MET/sqrt(HT)",				50, 0, 100);
	h_dphilw[i]   = new TH1D( Form( "dphilw_%s_%s"  , sampleName.Data(), regNames.at(i).Data()), "#Delta#phi (lep,W)",		50, 0, 3.5);
	h_njets[i]    = new TH1D( Form( "njets_%s_%s"   , sampleName.Data(), regNames.at(i).Data()), "Number of jets",            16, -0.5, 15.5);
	h_nbtags[i]   = new TH1D( Form( "nbtags_%s_%s"  , sampleName.Data(), regNames.at(i).Data()), "Number of b-tags",          7, -0.5, 6.5);
	h_ptj1[i]     = new TH1D( Form( "ptj1_%s_%s"    , sampleName.Data(), regNames.at(i).Data()), "Leading jet p_{T}",        40, 0, 1000);
	h_j1btag[i]   = new TH1D( Form( "j1btag_%s_%s"  , sampleName.Data(), regNames.at(i).Data()), "Is leading jet b-tagged?", 4, -0.5, 1.5);

	h_bgtype[i]->SetDirectory(rootdir);
	h_evttype[i]->SetDirectory(rootdir);
	h_sigyields[i]->SetDirectory(rootdir);

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
	h_ptj1[i]->SetDirectory(rootdir);
	h_j1btag[i]->SetDirectory(rootdir);

	TAxis* axis = h_bgtype[i]->GetXaxis();
	axis->SetBinLabel( 1, "ZtoNuNu" );
	axis->SetBinLabel( 2, "2+lep" );
	axis->SetBinLabel( 3, "1lepTop" );
	axis->SetBinLabel( 4, "1lepW" );
	axis->SetBinLabel( 5, "Other" );

	axis = h_evttype[i]->GetXaxis();
	axis->SetBinLabel( 1, "Data" );
	axis->SetBinLabel( 2, "Signals" );
	axis->SetBinLabel( 3, "ZtoNuNu" );
	axis->SetBinLabel( 4, "2+lep" );
	axis->SetBinLabel( 5, "1lepTop" );
	axis->SetBinLabel( 6, "1lepW" );

  }

  TH1D *h_yields = new TH1D( Form("srYields_%s", sampleName.Data()), "Yield by signal region", nSigRegs, 0.5, float(nSigRegs)+0.5);
  for( int i=0; i<nSigRegs; i++ ) h_yields->GetXaxis()->SetBinLabel( i+1, regNames.at(i) );
  h_yields->SetDirectory(rootdir);


  float yield_total = 0;
  float yield_unique = 0;
  float yield_filter = 0;
  float yield_vtx = 0;
  float yield_1goodlep = 0;
  // float yield_lepSel = 0;
  float yield_2lepveto = 0;
  float yield_trkVeto = 0;
  float yield_tauVeto = 0;
  float yield_njets = 0;
  float yield_1bjet = 0;
  float yield_METcut = 0;
  float yield_MTcut = 0;
  float yield_dPhi = 0;
  float yield_chi2 = 0;

  int yGen_total = 0;
  int yGen_unique = 0;
  int yGen_filter = 0;
  int yGen_vtx = 0;
  int yGen_1goodlep = 0;
  // int yGen_lepSel = 0;
  int yGen_2lepveto = 0;
  int yGen_trkVeto = 0;
  int yGen_tauVeto = 0;
  int yGen_njets = 0;
  int yGen_1bjet = 0;
  int yGen_METcut = 0;
  int yGen_MTcut = 0;
  int yGen_dPhi = 0;
  int yGen_chi2 = 0;

  ////////////////////////////////////////////////////////////////////
  // Bad event filters and such for data

  eventFilter badEventFilter;
  if( mySample->IsData() ) {
	cout << "Loading bad event files..." << endl;
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/eventlist_DoubleEG_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/eventlist_DoubleMuon_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/eventlist_HTMHT_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/eventlist_JetHT_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/eventlist_MET_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/eventlist_MuonEG_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/eventlist_SingleElectron_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/eventlist_SingleMuon_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/eventlist_SinglePhoton_csc2015.txt");
    // new lists: supposed to include above but do not always
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/DoubleEG_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/DoubleMuon_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/HTMHT_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/JetHT_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/MET_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/MuonEG_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/SingleElectron_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/SingleMuon_csc2015.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/SinglePhoton_csc2015.txt");
    // not all samples have events which failed the ecal SC filter
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/DoubleEG_ecalscn1043093.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/HTMHT_ecalscn1043093.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/JetHT_ecalscn1043093.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/MET_ecalscn1043093.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/SinglePhoton_ecalscn1043093.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/SingleElectron_ecalscn1043093.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/SingleMuon_ecalscn1043093.txt");
    // Some new filters pointed out by HJ on Feb-2
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/csc2015_Dec01.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/ecalscn1043093_Dec01.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/badResolutionTrack_Jan13.txt");
    badEventFilter.loadBadEventList("/nfs-6/userdata/mt2utils/muonBadTrack_Jan13.txt");
    cout << " ... finished!" << endl;
  }


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
	  if(      sampleName == "tt2l"  && genlepsfromtop() != 2 ) continue;  //Require 2 leps from top in "tt2l" events
	  else if( sampleName == "tt1l"  && genlepsfromtop() != 1 ) continue;  //Require 1 lep from top in "tt1l" events


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
		evtWeight = myAnalysis->GetLumi() * 1000. * xsec() / nEvtsSample;
	  }
	  else evtWeight = myAnalysis->GetLumi() * scale1fb();


	  // Count the number of events processed
	  yield_total += evtWeight;
	  yGen_total++;

	  // Remove duplicate events in data
	  if( is_data() ) {
		duplicate_removal::DorkyEventIdentifier id( run(), evt(), ls() );
        if( is_duplicate(id) ) continue;
		yield_unique += evtWeight;
		yGen_unique++;
      }

	  // MET filters and bad event filters for data
	  if( is_data() ) {
		if( badEventFilter.eventFails( run(), ls(), evt() ) ) continue;
		if( !filt_cscbeamhalo() ) continue;
		if(	!filt_eebadsc() ) continue;
		if( !filt_goodvtx() ) continue;
		if( !filt_hbhenoise() ) continue;
		yield_filter += evtWeight;
		yGen_filter++;
	  }

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

	  // Second lepton veto
	  if( nvetoleps() > 1 && ROOT::Math::VectorUtil::DeltaR( lep1_p4(), lep2_p4() ) > 0.01 ) continue;
	  yield_2lepveto += evtWeight;
	  yGen_2lepveto++;

	  // Track veto
	  if( !PassTrackVeto_v3() ) continue;
	  yield_trkVeto += evtWeight;
	  yGen_trkVeto++;

	  // Tau veto
	  if( !PassTauVeto() ) continue;
	  yield_tauVeto += evtWeight;
	  yGen_tauVeto++;

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

	  j1pt = ak4pfjets_pt().at(0);

	  // B-tag requirement
	  if( ngoodbtags() < 1 ) continue;
	  yield_1bjet += evtWeight;
	  yGen_1bjet++;

	  j1_isBtag = ak4pfjets_passMEDbtag().at(0);

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
	  if(   isZtoNuNu() )        category = 1;   // Z to nu nu
	  else if( is2lep() )        category = 2;   // 2 or more leptons
	  else if( is1lepFromTop() ) category = 3;   // 1 lepton from top quark
	  else if( is1lepFromW() )   category = 4;   // 1 lepton from a W not from top
	  else                       category = 5;   // Other

	  int evtType = -99;
	  if(      mySample->IsData()   ) evtType = 1;
	  else if( mySample->IsSignal() ) evtType = 2;
	  else                            evtType = 2+category;

	  ///////////////////////////////////////////
	  // Signal region cuts and histo filling

	  // If the event passes the SR cuts, store which background type this event is, and fill histograms
	  for( int i=0; i<nSigRegs; i++ ) {

		if( !sigRegions.at(i).PassAllCuts() ) continue;

		h_bgtype[i]->Fill( category,                    evtWeight );
		h_evttype[i]->Fill( evtType,                    evtWeight );

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
		h_ptj1[i]->Fill(    j1pt,                       evtWeight );
		h_j1btag[i]->Fill(  j1_isBtag,                  evtWeight );

		h_yields->Fill(     float(i+1),                 evtWeight );

		if( mySample->IsSignal() ) h_sigyields[i]->Fill( mass_stop(), mass_lsp(), evtWeight );
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
  if( mySample->IsData() ) {
	printf("Events passing duplicate removal:   %10.2f %9i\n", yield_unique   , yGen_unique       );
	printf("Events passing event filters:       %10.2f %9i\n", yield_filter   , yGen_filter       );
  }
  printf("Events with 1st vertex good:        %10.2f %9i\n", yield_vtx		, yGen_vtx			);
  printf("Events with at least 1 good lepton: %10.2f %9i\n", yield_1goodlep	, yGen_1goodlep		);
  // printf("Events passing lepton selection:    %10.2f %9i\n", yield_lepSel	, yGen_lepSel		);
  printf("Events passing second lepton veto:  %10.2f %9i\n", yield_2lepveto	, yGen_2lepveto		);
  printf("Events passing track veto:          %10.2f %9i\n", yield_trkVeto	, yGen_trkVeto		);
  printf("Events passing tau veto:            %10.2f %9i\n", yield_tauVeto	, yGen_tauVeto		);
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
	  if( h_evttype[j]->GetBinContent(k+2) < 0.0 ) {
		h_evttype[j]->SetBinContent(k+2, 0.);
		h_evttype[j]->SetBinError(k+2, 0.);
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
  TFile* plotfile = new TFile( myAnalysis->GetFileName() , "UPDATE");
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

	// Build up histo of signal yields
	if( mySample->IsSignal() ) {
	  TH2D* hTemp2 = (TH2D*)plotfile->Get( h_sigyields[j]->GetName() );
	  if( hTemp2 != 0 ) h_sigyields[j]->Add( hTemp2 );
	  h_sigyields[j]->Write( "", TObject::kOverwrite );
	}

	// Build up histo of yields by bkg type
	TH1D* hTemp = (TH1D*)plotfile->Get( h_evttype[j]->GetName() );
	if( hTemp != 0 ) h_evttype[j]->Add( hTemp );
	h_evttype[j]->Write( "", TObject::kOverwrite );

	h_ptj1[j]->Write();
	h_j1btag[j]->Write();
  }
  h_yields->Write();

  plotfile->Close();

  for( int j=0; j<nSigRegs; j++ ) {
	delete h_bgtype[j];
	delete h_evttype[j];
	delete h_sigyields[j];
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
	delete h_ptj1[j];
	delete h_j1btag[j];
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
