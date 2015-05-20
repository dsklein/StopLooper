{

  gROOT->ProcessLine(".L ScanChain.C+");

  cout << "ttbar 1l" << endl;
  TChain *ch_tt1l = new TChain("t"); 
  ch_tt1l->Add("../babies/ttbar.root");
  // ch_tt1l->Add("../babies/tbar_sch.root");
  // ch_tt1l->Add("../babies/tbar_tch.root");
  // ch_tt1l->Add("../babies/t_sch.root");
  // ch_tt1l->Add("../babies/t_tch.root");
  // ch_tt1l->Add("../babies/t_tW.root");
  // ch_tt1l->Add("../babies/tbar_tW.root");
  ScanChain(ch_tt1l, 1);

  cout << "ttbar 2l" << endl;
  TChain *ch_tt2l = new TChain("t"); 
  ch_tt2l->Add("../babies/ttbar.root");
  ScanChain(ch_tt2l, 2);

  cout << "Stop 850 100" << endl;
  TChain *ch_stop_850 = new TChain("t"); 
  ch_stop_850->Add("../babies/stop_850_100.root");
  ScanChain(ch_stop_850);
}
