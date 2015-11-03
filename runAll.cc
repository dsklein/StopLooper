#include "analysis.h"
#include "sample.h"
#include "runAll.h"

#include "TChain.h"
#include "TFile.h"
#include "TString.h"

int main() {



  // Clear out the file to be filled with histograms  
  TFile* outfile = new TFile("plots.root", "RECREATE");
  outfile->Close();

  TString bkgPath = "/nfs-7/userdata/stopRun2/StopBabies__CMS3_V07-04-XX/Spring15_25ns_Samples/StopBabyMaker__v7.4.x_v5/Skim__METge30__LEPge1_elPt20_elEta2p1_muPt20_muEta2p1_vetoElPt5_vetoElEta2p4_vetoMuPt5_vetoMuEta2p4__JETge1_jPt30_jEta2p4__20150728/";

  TString sigPath = "/nfs-7/userdata/stopRun2/StopBabies__CMS3_V07-04-XX/Phys14Signals/StopBabyMaker__v7.4.x_v5/Skim__METge30__LEPge1_elPt20_elEta2p1_muPt20_muEta2p1_vetoElPt5_vetoElEta2p4_vetoMuPt5_vetoMuEta2p4__JETge1_jPt30_jEta2p4__20150728/";

  // Signal samples

  TChain *ch_stop850 = new TChain("t");
  ch_stop850->Add( sigPath + "stop_850_100.root");

  TChain *ch_stop650 = new TChain("t");
  ch_stop650->Add( sigPath + "stop_650_325.root");

  TChain *ch_stop500 = new TChain("t");
  ch_stop500->Add( sigPath + "stop_500_325.root");

  TChain *ch_stop425 = new TChain("t");
  ch_stop425->Add( sigPath + "stop_425_325.root");


  // Background samples

  TChain *ch_ttbar = new TChain("t");
  ch_ttbar->Add( bkgPath + "ttbar_powheg_pythia8_25ns.root" );

  TChain *ch_wjets = new TChain("t");
  ch_wjets->Add( bkgPath + "WJetsToLNu_HT100To200_madgraph_pythia8_25ns.root" );
  ch_wjets->Add( bkgPath + "WJetsToLNu_HT200To400_madgraph_pythia8_25ns.root" );
  ch_wjets->Add( bkgPath + "WJetsToLNu_HT400To600_madgraph_pythia8_25ns.root" );
  ch_wjets->Add( bkgPath + "WJetsToLNu_HT600ToInf_madgraph_pythia8_25ns.root" );

  TChain *ch_dy = new TChain("t");
  ch_dy->Add( bkgPath + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns.root" );
  ch_dy->Add( bkgPath + "DYJetsToLL_m50_amcnlo_pythia8_25ns.root" );

  TChain *ch_stch = new TChain("t");
  ch_stch->Add( bkgPath + "t_sch_4f_amcnlo_pythia8_25ns.root" );
  ch_stch->Add( bkgPath + "t_tch_4f_powheg_pythia8_25ns.root" );
  ch_stch->Add( bkgPath + "tbar_tch_4f_powheg_pythia8_25ns.root" );

  TChain *ch_sttw = new TChain("t");
  ch_sttw->Add( bkgPath + "t_tW_5f_powheg_pythia8_25ns.root" );
  ch_sttw->Add( bkgPath + "t_tbarW_5f_powheg_pythia8_25ns.root" );

  TChain *ch_ttw = new TChain("t");
  ch_ttw->Add( bkgPath + "TTWJetsToLNu_amcnlo_pythia8_25ns.root" );

  TChain *ch_ttz = new TChain("t");
  ch_ttz->Add( bkgPath + "TTZToLLNuNu_M-10_amcnlo_pythia8_25ns.root" );
  ch_ttz->Add( bkgPath + "TTZToQQ_amcnlo_pythia8_25ns.root" );

  TChain *ch_tzq = new TChain("t");
  ch_tzq->Add( bkgPath + "tZq_ll_4f_amcnlo_pythia8_25ns.root" );
  ch_tzq->Add( bkgPath + "tZq_nunu_4f_amcnlo_pythia8_25ns.root" );

  TChain *ch_vv = new TChain("t");
  ch_vv->Add( bkgPath + "WWTo2l2Nu_powheg_25ns.root" );
  ch_vv->Add( bkgPath + "WWToLNuQQ_powheg_25ns.root" );
  ch_vv->Add( bkgPath + "WZTo3LNu_powheg_pythia8_25ns.root" );
  ch_vv->Add( bkgPath + "WZTo2L2Q_amcnlo_pythia8_25ns.root" );
  ch_vv->Add( bkgPath + "WZTo1Lnu2Q_amcnlo_pythia8_25ns.root" );
  ch_vv->Add( bkgPath + "ZZTo4L_powheg_pythia8_25ns.root" );
  ch_vv->Add( bkgPath + "ZZTo2L2Q_amcnlo_pythia8_25ns.root" );
  ch_vv->Add( bkgPath + "ZZTo2L2Nu_powheg_pythia8_25ns.root" );
  ch_vv->Add( bkgPath + "ZZTo2Q2Nu_amcnlo_pythia8_25ns.root" );


  TChain *ch_singletop = new TChain("t");
  ch_singletop->Add( ch_stch );
  ch_singletop->Add( ch_sttw );

  TChain *ch_rare = new TChain("t");
  ch_rare->Add( ch_ttw );
  ch_rare->Add( ch_ttz );
  ch_rare->Add( ch_vv );
  ch_rare->Add( ch_tzq );


  //////////////////////////////////////////////////////////////////////////
  // Run ScanChain
  ScanChain(ch_stop850, "stop850");
  ScanChain(ch_stop650, "stop650");
  ScanChain(ch_stop500, "stop500");
  ScanChain(ch_stop425, "stop425");  
  ScanChain(ch_ttbar, "tt2l");
  ScanChain(ch_ttbar, "tt1l"); //Same baby, pick out different final state
  // ScanChain(ch_stch, "STstchan");
  // ScanChain(ch_sttw, "STtWchan");
  ScanChain(ch_singletop, "singletop");
  // ScanChain(ch_wjets, "Wb");
  // ScanChain(ch_wjets, "Wucsd"); //Same baby, pick out different final state
  ScanChain(ch_wjets, "wjets"); //
  ScanChain(ch_dy, "dy");
  // ScanChain(ch_ttw, "ttw");
  // ScanChain(ch_ttz, "ttz");
  // ScanChain(ch_vv, "vv");
  ScanChain(ch_rare, "rare");


  ////////////////////////////////////////////////////////////////////////////
  // Make "analysis" object out of "samples", and pass them to makeTables and makeStack

  sample* stop850   = new sample( "stop850", "T2tt (850,100)", kBlue+3,   true );
  sample* stop650   = new sample( "stop650", "T2tt (650,325)", kGreen+3,  true );
  sample* stop500   = new sample( "stop500", "T2tt (500,325)", kMagenta+3,true );
  sample* stop425   = new sample( "stop425", "T2tt (425,325)", kOrange+7, true );

  sample* tt2l      = new sample( "tt2l", "$t\\bar{t} \\rightarrow 2l$", "t#bar{t} #rightarrow 2l", kGreen-4,  false );
  sample* tt1l      = new sample( "tt1l", "$t\\bar{t} \\rightarrow 1l$", "t#bar{t} #rightarrow 1l", kCyan-3,   false );
  sample* wjets     = new sample( "wjets",   "W+Jets",         kRed-7,    false );
  sample* dy        = new sample( "dy",      "Drell-Yan",      kOrange-2, false );
  sample* singletop = new sample( "singletop", "Single Top",   kMagenta-5,false );
  sample* rare      = new sample( "rare",    "Rare",           kRed+2,    false );

  analysis* ThisAnalysis = new analysis;
  ThisAnalysis->AddSample(stop850);
  ThisAnalysis->AddSample(stop650);
  ThisAnalysis->AddSample(stop500);
  ThisAnalysis->AddSample(stop425);
  ThisAnalysis->AddSample(tt2l);
  ThisAnalysis->AddSample(tt1l);
  ThisAnalysis->AddSample(wjets);
  ThisAnalysis->AddSample(dy);
  ThisAnalysis->AddSample(singletop);
  ThisAnalysis->AddSample(rare);


  makeTables( ThisAnalysis );
  makeStack( ThisAnalysis );

  return 0;
}
