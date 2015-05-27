{

  gROOT->ProcessLine(".L ScanChain.C+");
  gROOT->ProcessLine(".L dataMCplotMaker.cc+");
  gROOT->ProcessLine(".L makeStack.C+");

  TFile* outfile = new TFile("plots.root", "RECREATE");
  outfile->Close();

  cout << "single top" << endl;
  TChain *ch_t1l = new TChain("t"); 
  // ch_t1l->Add("../babies/ttbar.root");
  ch_t1l->Add("../babies/tbar_sch.root");
  ch_t1l->Add("../babies/tbar_tch.root");
  ch_t1l->Add("../babies/t_sch.root");
  ch_t1l->Add("../babies/t_tch.root");
  ch_t1l->Add("../babies/t_tW.root");
  ch_t1l->Add("../babies/tbar_tW.root");
  ScanChain(ch_t1l, "top1l");  // May or may not want to require exactly 1 gen lepton

  cout << "ttbar 2l" << endl;
  TChain *ch_tt2l = new TChain("t"); 
  ch_tt2l->Add("../babies/ttbar.root");
  ScanChain(ch_tt2l, "tt2l", 2);

  cout << "ttbar 1l" << endl;
  ScanChain(ch_tt2l ,"tt1l", 1);

  cout << "W+jets" << endl;
  TChain *ch_wjets = new TChain("t"); 
  ch_wjets->Add("../babies/wjets.root");
  ScanChain(ch_wjets, "wjets");

  cout << "Rare" << endl;
  TChain *ch_rare = new TChain("t");
  ch_rare->Add("../babies/ttwjets.root");
  ch_rare->Add("../babies/ttzjets.root");
  ch_rare->Add("../babies/wzjets.root");
  ch_rare->Add("../babies/zz.root");
  ch_rare->Add("../babies/ttbar.root");
  ScanChain(ch_rare, "rare", 0);


  cout << "Stop 425 325" << endl;
  TChain *ch_stop_425 = new TChain("t"); 
  ch_stop_425->Add("../babies/stop_425_325.root");
  ScanChain(ch_stop_425, "stop425");

  cout << "Stop 500 325" << endl;
  TChain *ch_stop_500 = new TChain("t"); 
  ch_stop_500->Add("../babies/stop_500_325.root");
  ScanChain(ch_stop_500, "stop500");

  cout << "Stop 650 325" << endl;
  TChain *ch_stop_650 = new TChain("t"); 
  ch_stop_650->Add("../babies/stop_650_325.root");
  ScanChain(ch_stop_650, "stop650");

  cout << "Stop 850 100" << endl;
  TChain *ch_stop_850 = new TChain("t"); 
  ch_stop_850->Add("../babies/stop_850_100.root");
  ScanChain(ch_stop_850, "stop850");

  makeStack();
}
