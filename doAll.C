{

  gROOT->ProcessLine(".L ScanChain.C+");
  gROOT->ProcessLine(".L makeTables.C+");
  gROOT->ProcessLine(".L dataMCplotMaker.cc+");
  gROOT->ProcessLine(".L makeStack.C+");

  // Clear out the file to be filled with histograms  
  TFile* outfile = new TFile("plots.root", "RECREATE");
  outfile->Close();

  TString babyPath = "/nfs-7/userdata/stopRun2/StopBabies__CMS3_V07-04-XX/Spring15_25ns_Samples/StopBabyMaker__v7.4.x_v5/Skim__METge30__LEPge1_elPt20_elEta2p1_muPt20_muEta2p1_vetoElPt5_vetoElEta2p4_vetoMuPt5_vetoMuEta2p4__JETge1_jPt30_jEta2p4__20150728/";



  TChain *ch_ttbar = new TChain("t");
  ch_ttbar->Add( babyPath + "ttbar_powheg_pythia8_25ns.root" );
  ScanChain(ch_ttbar, "tt2l");

  ScanChain(ch_ttbar, "tt1l"); //Same baby, pick out different final state

  TChain *ch_wjets = new TChain("t");
  ch_wjets->Add( babyPath + "WJetsToLNu_HT100To200_madgraph_pythia8_25ns.root" );
  ch_wjets->Add( babyPath + "WJetsToLNu_HT200To400_madgraph_pythia8_25ns.root" );
  ch_wjets->Add( babyPath + "WJetsToLNu_HT400To600_madgraph_pythia8_25ns.root" );
  ch_wjets->Add( babyPath + "WJetsToLNu_HT600ToInf_madgraph_pythia8_25ns.root" );
  ScanChain(ch_wjets, "Wb");

  ScanChain(ch_wjets, "Wucsd"); //Same baby, pick out different final state

  TChain *ch_dy = new TChain("t");
  ch_dy->Add( babyPath + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns.root" );
  ch_dy->Add( babyPath + "DYJetsToLL_m50_amcnlo_pythia8_25ns.root" );
  ScanChain(ch_dy, "dy");

  TChain *ch_ttw = new TChain("t");
  ch_ttw->Add( babyPath + "TTWJetsToLNu_amcnlo_pythia8_25ns.root" );
  ScanChain(ch_ttw, "ttw", 1);

  TChain *ch_ttz = new TChain("t");
  ch_ttz->Add( babyPath + "TTZToLLNuNu_M-10_amcnlo_pythia8_25ns.root" );
  ch_ttz->Add( babyPath + "TTZToQQ_amcnlo_pythia8_25ns.root" );
  ScanChain(ch_ttz, "ttz", 1);

  TChain *ch_stch = new TChain("t");
  // ch_stch->Add( babyPath + "t_sch_4f_amcnlo_pythia8_25ns.root" ); // Not a default sample
  ch_stch->Add( babyPath + "t_tch_4f_powheg_pythia8_25ns.root" );
  ch_stch->Add( babyPath + "tbar_tch_4f_powheg_pythia8_25ns.root" );
  ScanChain(ch_stch, "STstchan");

  TChain *ch_sttw = new TChain("t");
  ch_sttw->Add( babyPath + "t_tW_5f_powheg_pythia8_25ns.root" );
  ch_sttw->Add( babyPath + "t_tbarW_5f_powheg_pythia8_25ns.root" );
  ScanChain(ch_sttw, "STtWchan");

  TChain *ch_vv = new TChain("t");
  ch_vv->Add( babyPath + "WWTo2l2Nu_powheg_25ns.root" );
  // ch_vv->Add( babyPath + "WWToLNuQQ_powheg_25ns.root" ); // Not a default sample
  ch_vv->Add( babyPath + "WZ_pythia8_25ns.root" );
  ch_vv->Add( babyPath + "ZZ_pythia8_25ns.root" );
  ScanChain(ch_vv, "vv");


  makeTables();
  makeStack();
}
