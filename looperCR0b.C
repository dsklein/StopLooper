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
#include "../StopAnalysis_74x/CORE/Tools/badEventFilter.h"
#include "../StopAnalysis_74x/CORE/Tools/dorky/dorky.h"

// Custom
#include "analysis.h"
#include "sample.h"

using namespace std;
using namespace tas;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

// Global variables used in defining signal regions
extern bool j1_isBtag;
extern double j1pt;


int looperCR0b( analysis* myAnalysis, sample* mySample, int nEvents = -1, bool fast = true) {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Setup
  TChain *chain = mySample->GetChain();
  TString sampleName = mySample->GetLabel();
  const int nSigRegs = myAnalysis->GetSigRegionsAll().size();
  const int nVariations = myAnalysis->GetSystematics(false).size();
  bool isFastsim = mySample->IsSignal();
  cout << "\nSample: " << sampleName.Data() << endl;

  /////////////////////////////////////////////////////////
  // Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  TH1::SetDefaultSumw2();

  TH1D* h_bkgtype[nSigRegs][nVariations+1];
  TH1D* h_evttype[nSigRegs][nVariations+1];
  TH2D* h_sigyields[nSigRegs][nVariations+1];

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
  TH1D *h_modtop[nSigRegs];

  vector<TString> regNames = myAnalysis->GetSigRegionLabelsAll();
  vector<sigRegion> sigRegions = myAnalysis->GetSigRegionsAll();
  vector<systematic*> variations = myAnalysis->GetSystematics(false);

  for( int i=0; i<nSigRegs; i++ ) {

	TString plotLabel = sampleName + "_" + regNames.at(i);

	h_bkgtype[i][0]  = new TH1D(  "bkgtype_" + plotLabel, "Yield by background type",  5, 0.5, 5.5);
	h_evttype[i][0]  = new TH1D(  "evttype_" + regNames.at(i), "Yield by event type",  6, 0.5, 6.5);
	h_sigyields[i][0] = new TH2D( "sigyields_" + regNames.at(i), "Signal yields by mass point", 37, 87.5, 1012.5, 21, -12.5, 512.5 );

	for( int j=1; j<=nVariations; j++ ) {
	  TString varName = variations.at(j-1)->GetNameLong();
	  h_bkgtype[i][j]   = new TH1D( "bkgtype_" + plotLabel + "_" + varName, "Yield by background type",  5, 0.5, 5.5);
	  h_evttype[i][j]   = new TH1D( "evttype_" + regNames.at(i) + "_" + varName, "Yield by event type",  6, 0.5, 6.5);
	  h_sigyields[i][j] = new TH2D( "sigyields_" + regNames.at(i) + "_" + varName, "Signal yields by mass point", 37, 87.5, 1012.5, 21, -12.5, 512.5 );
	}

	h_mt[i]       = new TH1D(  "mt_"      + plotLabel, "Transverse mass",			80, 0, 800);
	h_met[i]      = new TH1D(  "met_"     + plotLabel, "MET",						40, 0, 1000);
	h_mt2w[i]     = new TH1D(  "mt2w_"    + plotLabel, "MT2W",						50, 0, 500);
	h_chi2[i]     = new TH1D(  "chi2_"    + plotLabel, "Hadronic #chi^{2}", 		50, 0, 15);
	h_htratio[i]  = new TH1D(  "htratio_" + plotLabel, "H_{T} ratio",				50, 0, 1);
	h_mindphi[i]  = new TH1D(  "mindphi_" + plotLabel, "min #Delta#phi(j12,MET)",	50, 0, 4);
	h_ptb1[i]     = new TH1D(  "ptb1_"    + plotLabel, "p_{T} (b1)",				100, 0, 500);
	h_drlb1[i]    = new TH1D(  "drlb1_"   + plotLabel, "#DeltaR (lep, b1)", 		50, 0, 5);
	h_ptlep[i]    = new TH1D(  "ptlep_"   + plotLabel, "p_{T} (lep)",				100, 0, 500);
	h_metht[i]    = new TH1D(  "metht_"   + plotLabel, "MET/sqrt(HT)",				50, 0, 100);
	h_dphilw[i]   = new TH1D(  "dphilw_"  + plotLabel, "#Delta#phi (lep,W)",		50, 0, 3.5);
	h_njets[i]    = new TH1D(  "njets_"   + plotLabel, "Number of jets",            16, -0.5, 15.5);
	h_nbtags[i]   = new TH1D(  "nbtags_"  + plotLabel, "Number of b-tags",          7, -0.5, 6.5);
	h_ptj1[i]     = new TH1D(  "ptj1_"    + plotLabel, "Leading jet p_{T}",        40, 0, 1000);
	h_j1btag[i]   = new TH1D(  "j1btag_"  + plotLabel, "Is leading jet b-tagged?", 2, -0.5, 1.5);
	h_modtop[i]   = new TH1D(  "modtop_"  + plotLabel, "Modified topness",         30, -15., 15.);

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
	h_modtop[i]->SetDirectory(rootdir);

	for( int j=0; j<=nVariations; j++ ) {
	  h_bkgtype[i][j]->SetDirectory(rootdir);
	  h_evttype[i][j]->SetDirectory(rootdir);
	  h_sigyields[i][j]->SetDirectory(rootdir);

	  TAxis* axis = h_bkgtype[i][j]->GetXaxis();
	  axis->SetBinLabel( 1, "ZtoNuNu" );
	  axis->SetBinLabel( 2, "2+lep" );
	  axis->SetBinLabel( 3, "1lepTop" );
	  axis->SetBinLabel( 4, "1lepW" );
	  axis->SetBinLabel( 5, "Other" );

	  axis = h_evttype[i][j]->GetXaxis();
	  axis->SetBinLabel( 1, "Data" );
	  axis->SetBinLabel( 2, "Signals" );
	  axis->SetBinLabel( 3, "ZtoNuNu" );
	  axis->SetBinLabel( 4, "2+lep" );
	  axis->SetBinLabel( 5, "1lepTop" );
	  axis->SetBinLabel( 6, "1lepW" );
	}

  }

  TH1D *h_yields = new TH1D( Form("srYields_%s", sampleName.Data()), "Yield by signal region", nSigRegs, 0.5, float(nSigRegs)+0.5);
  for( int i=0; i<nSigRegs; i++ ) h_yields->GetXaxis()->SetBinLabel( i+1, regNames.at(i) );
  h_yields->SetDirectory(rootdir);


  double yield_total = 0;
  double yield_unique = 0;
  double yield_filter = 0;
  double yield_vtx = 0;
  double yield_1goodlep = 0;
  double yield_lepSel = 0;
  double yield_2lepveto = 0;
  double yield_trkVeto = 0;
  double yield_tauVeto = 0;
  double yield_njets = 0;
  double yield_0bjet = 0;
  double yield_METcut = 0;
  double yield_MTcut = 0;
  double yield_dPhi = 0;
  double yield_chi2 = 0;

  int yGen_total = 0;
  int yGen_unique = 0;
  int yGen_filter = 0;
  int yGen_vtx = 0;
  int yGen_1goodlep = 0;
  int yGen_lepSel = 0;
  int yGen_2lepveto = 0;
  int yGen_trkVeto = 0;
  int yGen_tauVeto = 0;
  int yGen_njets = 0;
  int yGen_0bjet = 0;
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

	duplicate_removal::clear_list();
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

	// Load event weight histograms
	TH2F* hNEvts = (TH2F*)file.Get("histNEvts");
	TH3D* hCounterSMS = (TH3D*)file.Get("h_counterSMS");
	TH1D* hCounter = (TH1D*)file.Get("h_counter");
	myHelper.Setup( isFastsim, hCounter, hNEvts, hCounterSMS );

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
	  if(      sampleName == "tt2l"  && gen_nfromtleps_() != 2 ) continue;  //Require 2 leps from top in "tt2l" events
	  else if( sampleName == "tt1l"  && gen_nfromtleps_() != 1 ) continue;  //Require 1 lep from top in "tt1l" events

	  //FastSim anomalous event filter
	  if( isFastsim && filt_fastsimjets() ) continue;


	  /////////////////////////////////
	  // Set event weight

	  double evtWeight = 1.;
	  double lepNorm = 1.;
	  double lepNorm_veto = 1.;
	  double lepNorm_FS = 1.;
	  double btagNorm = 1.;
	  double isrNorm = 1.;

	  if( is_data() || mySample->IsData() ) evtWeight = 1.;
	  else if( mySample->IsSignal() ) {
		double nEvtsSample = hNEvts->GetBinContent( hNEvts->FindBin( mass_stop(), mass_lsp() ) );
		int binx = hCounterSMS->GetXaxis()->FindBin( mass_stop() );
		int biny = hCounterSMS->GetYaxis()->FindBin( mass_lsp()  );
		lepNorm = nEvtsSample / hCounterSMS->GetBinContent(binx,biny,27);
		lepNorm_veto = nEvtsSample / hCounterSMS->GetBinContent(binx,biny,30);
		lepNorm_FS = nEvtsSample / hCounterSMS->GetBinContent(binx,biny,33);
		btagNorm = nEvtsSample / hCounterSMS->GetBinContent(binx,biny,14);
		isrNorm = nEvtsSample / hCounterSMS->GetBinContent(binx, biny, 19);
		evtWeight = myAnalysis->GetLumi() * 1000. * xsec() / nEvtsSample;
	  }
	  else {
		evtWeight = myAnalysis->GetLumi() * scale1fb();
		double nEvtsSample = hCounter->GetBinContent(22);
		lepNorm = nEvtsSample / hCounter->GetBinContent(28);
		lepNorm_veto = nEvtsSample / hCounter->GetBinContent(31);
		lepNorm_FS = nEvtsSample / hCounter->GetBinContent(34);
		btagNorm = nEvtsSample / hCounter->GetBinContent(14);
		isrNorm = nEvtsSample / hCounter->GetBinContent(19);
	  }

	  if( !is_data() ) {
	  	evtWeight *= weight_lepSF()     * lepNorm;
	  	evtWeight *= weight_vetoLepSF() * lepNorm_veto;
	  	evtWeight *= weight_btagsf() * btagNorm;
	  	if( isFastsim ) evtWeight *= weight_lepSF_fastSim() * lepNorm_FS;
		if( mySample->IsSignal()  ) evtWeight *= weight_ISR() * isrNorm;
	  }


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

	  // MET filters, bad event filters, and triggers for data
	  if( is_data() ) {
		if( badEventFilter.eventFails( run(), ls(), evt() ) ) continue;
		if( !filt_met() ) continue;
		if( !HLT_SingleEl() && !HLT_SingleMu() && !HLT_MET() ) continue;
		yield_filter += evtWeight;
		yGen_filter++;
	  }

	  // First vertex must be good
	  // if( firstGoodVtxIdx() != 0 ) continue;
	  yield_vtx += evtWeight;
	  yGen_vtx++;

	  // Must have at least 1 good lepton
	  if( ngoodleps() < 1 ) continue;
	  yield_1goodlep += evtWeight;
	  yGen_1goodlep++;

	  // Lep 1 must pass lepton selections
	  if( abs(lep1_pdgid())==11 ) {
	  	if( lep1_p4().pt() < 20. ) continue;
	  	if( fabs(lep1_p4().eta()) > 1.4442 ) continue;
	  	if( !lep1_passMediumID() ) continue;
	  }
	  else if( abs(lep1_pdgid())==13 ) {
	  	if( lep1_p4().pt() < 20. ) continue;
	  	if( fabs(lep1_p4().eta()) > 2.4 ) continue;
	  	if( !lep1_passTightID() ) continue;
	  }
	  yield_lepSel += evtWeight;
	  yGen_lepSel++;

	  // Second lepton veto
	  if( nvetoleps() > 1 && ROOT::Math::VectorUtil::DeltaR( lep1_p4(), lep2_p4() ) > 0.01 ) continue;
	  yield_2lepveto += evtWeight;
	  yGen_2lepveto++;

	  // Track veto
	  if( !PassTrackVeto() ) continue;
	  yield_trkVeto += evtWeight;
	  yGen_trkVeto++;

	  // Tau veto
	  if( !PassTauVeto() ) continue;
	  yield_tauVeto += evtWeight;
	  yGen_tauVeto++;

	  // N-jet requirement
	  if( ngoodjets() < 2 ) continue;
	  yield_njets += evtWeight;
	  yGen_njets++;

	  j1pt = ak4pfjets_p4().at(0).pt();

	  // B-tag requirement
	  if( ngoodbtags() != 0 ) continue;
	  yield_0bjet += evtWeight;
	  yGen_0bjet++;

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

	  // Quickly calculate some variables
	  double metSqHT = pfmet() / sqrt( ak4_HT() );

	  const TVector3 lepVec( lep1_p4().x(), lep1_p4().y(), lep1_p4().z() );
	  const TVector3 metVec( pfmet()*cos(pfmet_phi()), pfmet()*sin(pfmet_phi()), 0 );
	  const TVector3 wVec = lepVec + metVec;
	  double dPhiLepW = fabs( lepVec.DeltaPhi(wVec) );

	  double drLepLeadb = ROOT::Math::VectorUtil::DeltaR( lep1_p4(), ak4pfjets_leadMEDbjet_p4() );

	  ///////////////////////////////////////////
	  // Signal region cuts and histo filling

	  // If the event passes the SR cuts, store which background type this event is, and fill histograms
	  for( int i=0; i<nSigRegs; i++ ) {

		if( !sigRegions.at(i).PassAllCuts() ) continue;

		h_bkgtype[i][0]->Fill( category,                   evtWeight );
		h_evttype[i][0]->Fill( evtType,                    evtWeight );
		if( mySample->IsSignal() ) h_sigyields[i][0]->Fill( mass_stop(), mass_lsp(), evtWeight );

		h_mt[i]->Fill(      mt_met_lep(), 				evtWeight );
		h_met[i]->Fill(     pfmet(),					evtWeight );
		h_mt2w[i]->Fill(	MT2W(),   					evtWeight );
		h_chi2[i]->Fill(	hadronic_top_chi2(),		evtWeight );
		h_htratio[i]->Fill( ak4_htratiom(),				evtWeight );
		h_mindphi[i]->Fill( mindphi_met_j1_j2(),		evtWeight );
		h_ptb1[i]->Fill(	ak4pfjets_leadMEDbjet_p4().pt(),	evtWeight );
		h_drlb1[i]->Fill(   drLepLeadb,					evtWeight );
		h_ptlep[i]->Fill(   lep1_p4().pt(),				evtWeight );
		h_metht[i]->Fill(   metSqHT,					evtWeight );
		h_dphilw[i]->Fill(  dPhiLepW,					evtWeight );
		h_njets[i]->Fill(   ngoodjets(),                evtWeight );
		h_nbtags[i]->Fill(  ngoodbtags(),               evtWeight );
		h_ptj1[i]->Fill(    j1pt,                       evtWeight );
		h_j1btag[i]->Fill(  j1_isBtag,                  evtWeight );
		h_modtop[i]->Fill(  topnessMod(),               evtWeight );

		h_yields->Fill(     double(i+1),                evtWeight );

		// Special systematic variation histograms
		for( int j=1; j<=nVariations; j++ ) {
		  h_bkgtype[i][j]->Fill( category,   evtWeight * variations.at(j-1)->GetWeight() );
		  h_evttype[i][j]->Fill( evtType,    evtWeight * variations.at(j-1)->GetWeight() );
		  if( mySample->IsSignal() ) h_sigyields[i][j]->Fill( mass_stop(), mass_lsp(), evtWeight * variations.at(j-1)->GetWeight() );
		}

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
	printf("Events passing filters and trigger: %10.2f %9i\n", yield_filter   , yGen_filter       );
  }
  printf("Events with 1st vertex good:        %10.2f %9i\n", yield_vtx		, yGen_vtx			);
  printf("Events with at least 1 good lepton: %10.2f %9i\n", yield_1goodlep	, yGen_1goodlep		);
  printf("Events passing lepton selection:    %10.2f %9i\n", yield_lepSel	, yGen_lepSel		);
  printf("Events passing second lepton veto:  %10.2f %9i\n", yield_2lepveto	, yGen_2lepveto		);
  printf("Events passing track veto:          %10.2f %9i\n", yield_trkVeto	, yGen_trkVeto		);
  printf("Events passing tau veto:            %10.2f %9i\n", yield_tauVeto	, yGen_tauVeto		);
  printf("Events with at least 2 jets:        %10.2f %9i\n", yield_njets	, yGen_njets		);
  printf("Events with zero b-tags:            %10.2f %9i\n", yield_0bjet	, yGen_0bjet		);
  printf("Events with MET > 250 GeV:          %10.2f %9i\n", yield_METcut	, yGen_METcut		);
  printf("Events with MT > 150 GeV:           %10.2f %9i\n", yield_MTcut	, yGen_MTcut		);
  printf("Events with min dPhi > 0.8:         %10.2f %9i\n", yield_dPhi		, yGen_dPhi			);
  // printf("Events with chi2 < 10:              %10.2f %9i\n", yield_chi2 	, yGen_chi2 		);
  printf("Yield after preselection:           %10.2f %9i\n", yield_chi2		, yGen_chi2			);

  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  // Zero negative values in each signal region
  for( int i=0; i<nSigRegs; i++ ) {
	for( int j=0; j<=nVariations; j++ ) {
	  bool negsFound = false;

	  // First zero any decay modes with negative yields
	  for( int k=1; k<= h_bkgtype[i][j]->GetNbinsX(); k++ ) {
		if( h_bkgtype[i][j]->GetBinContent(k) < 0.0 ) {
		  h_bkgtype[i][j]->SetBinContent(k, 0.);
		  h_bkgtype[i][j]->SetBinError(k, 0.);
		  negsFound = true;
		}
		if( h_evttype[i][j]->GetBinContent(k+2) < 0.0 ) {
		  h_evttype[i][j]->SetBinContent(k+2, 0.);
		  h_evttype[i][j]->SetBinError(k+2, 0.);
		}
	  }
	  // If any negative yields were found in any decay mode, recalculate the total yield
	  if( j==0 && negsFound ) {
		double newYield, newErr;
		newYield = h_bkgtype[i][0]->IntegralAndError( 0, -1, newErr );
		h_yields->SetBinContent(i+1, newYield);
		h_yields->SetBinError(i+1, newErr);
	  }
	}
  }

  // Store histograms and clean them up
  TFile* plotfile = new TFile( myAnalysis->GetPlotFileName(), "UPDATE");
  plotfile->cd();

  for( int i=0; i<nSigRegs; i++ ) {
	h_bkgtype[i][0]->Write();
	h_mt[i]->Write();
	h_met[i]->Write();
	h_mt2w[i]->Write();
	h_chi2[i]->Write();
	h_htratio[i]->Write();
	h_mindphi[i]->Write();
	h_ptb1[i]->Write();
	h_drlb1[i]->Write();
	h_ptlep[i]->Write();
	h_metht[i]->Write();
	h_dphilw[i]->Write();
	h_njets[i]->Write();
	h_nbtags[i]->Write();
	h_ptj1[i]->Write();
	h_j1btag[i]->Write();
	h_modtop[i]->Write();

	// Build up histo of signal yields
	if( mySample->IsSignal() ) {
	  TH2D* hTemp2 = (TH2D*)plotfile->Get( h_sigyields[i][0]->GetName() );
	  if( hTemp2 != 0 ) h_sigyields[i][0]->Add( hTemp2 );
	  h_sigyields[i][0]->Write( "", TObject::kOverwrite );
	}

	// Build up histo of yields by bkg type
	TH1D* hTemp = (TH1D*)plotfile->Get( h_evttype[i][0]->GetName() );
	if( hTemp != 0 ) h_evttype[i][0]->Add( hTemp );
	h_evttype[i][0]->Write( "", TObject::kOverwrite );
  }
  h_yields->Write();

  plotfile->Close();

  // Do similarly for the systematic variation histograms, but put them in a different file
  if( nVariations > 0 ) {
	TFile* systFile = new TFile(myAnalysis->GetSystFileName(), "UPDATE");
	systFile->cd();
	for( int i=0; i< nSigRegs; i++ ) {
	  for( int j=1; j<=nVariations; j++ ) {
		h_bkgtype[i][j]->Write();

		if( mySample->IsSignal() ) {
		  TH2D* hTemp2 = (TH2D*)systFile->Get( h_sigyields[i][j]->GetName() );
		  if( hTemp2 != 0 ) h_sigyields[i][j]->Add( hTemp2 );
		  h_sigyields[i][j]->Write( "", TObject::kOverwrite );
		}

		TH1D* hTemp = (TH1D*)systFile->Get( h_evttype[i][j]->GetName() );
		if( hTemp != 0 ) h_evttype[i][j]->Add( hTemp );
		h_evttype[i][j]->Write( "", TObject::kOverwrite );
	  }
	}
	systFile->Close();
  }

  // Cleanup
  for( int i=0; i<nSigRegs; i++ ) {
	for( int j=0; j<=nVariations; j++ ) {
	  delete h_bkgtype[i][j];
	  delete h_evttype[i][j];
	  delete h_sigyields[i][j];
	}
	delete h_mt[i];
	delete h_met[i];
	delete h_mt2w[i];
	delete h_chi2[i];
	delete h_htratio[i];
	delete h_mindphi[i];
	delete h_ptb1[i];
	delete h_drlb1[i];
	delete h_ptlep[i];
	delete h_metht[i];
	delete h_dphilw[i];
	delete h_njets[i];
	delete h_nbtags[i];
	delete h_ptj1[i];
	delete h_j1btag[i];
	delete h_modtop[i];
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
