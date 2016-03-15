#include "analysis.h"
#include "sample.h"
#include "runLooper.h"

#include "TChain.h"
#include "TFile.h"
#include "TString.h"

int main( int argc, char* argv[] ) {


  TString bkgPath = "/hadoop/cms/store/user/jgwood/condor/stop_1l_babies/stop_babies__CMS3_V070411__BabyMaker_V0704X_v9__20160127/merged_files/Skims_SR__20160202/";

  TString sigPath = "/nfs-7/userdata/stopRun2/signalbabies/";



  // Parse command-line argument(s)
  TString argument;
  if( argc==1 )       argument = "all";
  else if( argc > 1 ) argument = TString(argv[1]);

  if( argc > 2 || argument=="help" || argument=="h" ) {
	std::cout << "Usage: ./runLooper [arg]\n" << std::endl;

	std::cout << "Takes at most one argument from the following list:" << std::endl;
	std::cout << "help        show this message" << std::endl;
	std::cout << "[blank]     equivalent to 'all'" << std::endl;
	std::cout << "all         run ScanChain, makeTables, and makeStack" << std::endl;
	std::cout << "scanchain   run ScanChain only" << std::endl;
	std::cout << "plots       run makeStack only" << std::endl;
	std::cout << "tables      run makeTables only" << std::endl;
	std::cout << "output      run makeStack and makeTables only" << std::endl;
	return 0;
  }


  // Make Chains and run ScanChain
  if( argument=="all" || argument=="scanchain" || argument=="scan" || argument=="loop" ) {

	// Signal samples

	TChain *ch_stop700 = new TChain("t");
	ch_stop700->Add( sigPath + "Signal_T2tt_700_50.root");

	TChain *ch_stop600 = new TChain("t");
	ch_stop600->Add( sigPath + "Signal_T2tt_600_250.root");

	TChain *ch_stop300 = new TChain("t");
	ch_stop300->Add( sigPath + "Signal_T2tt_300_200.root");

	TChain *ch_stop275 = new TChain("t");
	ch_stop275->Add( sigPath + "Signal_T2tt_275_100.root");


	// Background samples

	TChain *ch_ttbar = new TChain("t");
	ch_ttbar->Add( bkgPath + "ttbar_powheg_pythia8_25ns_skimmed.root" );
	// ch_ttbar->Add( bkgPath + "ttbar_powheg_pythia8_ext3_25ns_skimmed.root" );

	TChain *ch_wjets = new TChain("t");
	ch_wjets->Add( bkgPath + "WJetsToLNu_HT100To200_madgraph_pythia8_25ns_skimmed.root" );
	ch_wjets->Add( bkgPath + "WJetsToLNu_HT200To400_madgraph_pythia8_25ns_skimmed.root" );
	ch_wjets->Add( bkgPath + "WJetsToLNu_HT400To600_madgraph_pythia8_25ns_skimmed.root" );
	ch_wjets->Add( bkgPath + "WJetsToLNu_HT600ToInf_madgraph_pythia8_25ns_skimmed.root" );

	TChain *ch_dy = new TChain("t");
	ch_dy->Add( bkgPath + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns_skimmed.root" );
	ch_dy->Add( bkgPath + "DYJetsToLL_m50_amcnlo_pythia8_25ns_skimmed.root" );

	TChain *ch_stch = new TChain("t");
	ch_stch->Add( bkgPath + "t_sch_4f_amcnlo_pythia8_25ns_skimmed.root" );
	ch_stch->Add( bkgPath + "t_tch_4f_powheg_pythia8_25ns_skimmed.root" );
	ch_stch->Add( bkgPath + "tbar_tch_4f_powheg_pythia8_25ns_skimmed.root" );

	TChain *ch_sttw = new TChain("t");
	ch_sttw->Add( bkgPath + "t_tW_5f_powheg_pythia8_25ns_skimmed.root" );
	ch_sttw->Add( bkgPath + "t_tbarW_5f_powheg_pythia8_25ns_skimmed.root" );

	TChain *ch_ttw = new TChain("t");
	ch_ttw->Add( bkgPath + "TTWJetsToLNu_amcnlo_pythia8_25ns_skimmed.root" );
	ch_ttw->Add( bkgPath + "TTWJetsToQQ_amcnlo_pythia8_25ns_skimmed.root" );

	TChain *ch_ttz = new TChain("t");
	ch_ttz->Add( bkgPath + "TTZToLLNuNu_M-10_amcnlo_pythia8_25ns_skimmed.root" );
	ch_ttz->Add( bkgPath + "TTZToQQ_amcnlo_pythia8_25ns_skimmed.root" );

	TChain *ch_tzq = new TChain("t");
	ch_tzq->Add( bkgPath + "tZq_ll_4f_amcnlo_pythia8_25ns_skimmed.root" );
	ch_tzq->Add( bkgPath + "tZq_nunu_4f_amcnlo_pythia8_25ns_skimmed.root" );

	TChain *ch_vv = new TChain("t");
	ch_vv->Add( bkgPath + "WWTo2l2Nu_powheg_25ns_skimmed.root" );
	ch_vv->Add( bkgPath + "WWToLNuQQ_powheg_25ns_skimmed.root" );
	ch_vv->Add( bkgPath + "WZTo3LNu_powheg_pythia8_25ns_skimmed.root" );
	ch_vv->Add( bkgPath + "WZTo2L2Q_amcnlo_pythia8_25ns_skimmed.root" );
	ch_vv->Add( bkgPath + "WZTo1LNu2Q_amcnlo_pythia8_25ns_skimmed.root" );
	ch_vv->Add( bkgPath + "ZZTo4L_powheg_pythia8_25ns_skimmed.root" );
	ch_vv->Add( bkgPath + "ZZTo2L2Q_amcnlo_pythia8_25ns_skimmed.root" );
	ch_vv->Add( bkgPath + "ZZTo2L2Nu_powheg_pythia8_25ns_skimmed.root" );
	ch_vv->Add( bkgPath + "ZZTo2Q2Nu_amcnlo_pythia8_25ns_skimmed.root" );


	TChain *ch_singletop = new TChain("t");
	ch_singletop->Add( ch_stch );
	ch_singletop->Add( ch_sttw );

	TChain *ch_rare = new TChain("t");
	ch_rare->Add( ch_ttw );
	ch_rare->Add( ch_ttz );
	ch_rare->Add( ch_vv );
	ch_rare->Add( ch_tzq );


	//////////////////////////////////////////////////////////////////////////
	// Reset output file and run ScanChain

	TFile* outfile = new TFile("plots.root", "RECREATE");
	outfile->Close();

	ScanChain(ch_stop700, "stop700");
	ScanChain(ch_stop600, "stop600");
	ScanChain(ch_stop300, "stop300");
	ScanChain(ch_stop275, "stop275");
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
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // Make "analysis" object out of "samples", and pass them to makeTables and makeStack

  //sample( "Storage name", "Display name", TColor, isData, isSignal

  sample* stop700   = new sample( "stop700", "T2tt (700,50)",  kBlue+3,   false, true );
  sample* stop600   = new sample( "stop600", "T2tt (600,250)", kGreen+3,  false, true );
  sample* stop300   = new sample( "stop300", "T2tt (300,200)", kMagenta+3,false, true );
  sample* stop275   = new sample( "stop275", "T2tt (275,100)", kOrange+7, false, true );

  sample* tt2l      = new sample( "tt2l", "$t\\bar{t} \\rightarrow 2l$", "t#bar{t} #rightarrow 2l", kCyan-3,   false, false );
  sample* tt1l      = new sample( "tt1l", "$t\\bar{t} \\rightarrow 1l$", "t#bar{t} #rightarrow 1l", kRed-7,    false, false );
  sample* singletop = new sample( "singletop", "Single Top",   kGreen-4,  false, false );
  sample* wjets     = new sample( "wjets",   "W+Jets",         kOrange-2, false, false );
  sample* dy        = new sample( "dy",      "Drell-Yan",      kRed+2,    false, false );
  sample* rare      = new sample( "rare",    "Rare",           kMagenta-5,false, false );

  analysis* ThisAnalysis = new analysis;
  ThisAnalysis->AddSample(stop700);
  ThisAnalysis->AddSample(stop600);
  ThisAnalysis->AddSample(stop300);
  ThisAnalysis->AddSample(stop275);
  ThisAnalysis->AddSample(tt2l);
  ThisAnalysis->AddSample(tt1l);
  ThisAnalysis->AddSample(singletop);
  ThisAnalysis->AddSample(wjets);
  ThisAnalysis->AddSample(dy);
  ThisAnalysis->AddSample(rare);

  std::vector<TString> compressed  = {"compr250",  "compr350"};
  std::vector<TString> boosted     = {"boost250",  "boost350"};
  std::vector<TString> lowDMreg    = {"low250",  "low325"};
  std::vector<TString> highDMreg = {"high250", "high350", "high450"};
  std::vector<TString> inclusive = {"inclusive"};
  ThisAnalysis->AddSigRegs( compressed );
  ThisAnalysis->AddSigRegs( boosted );
  ThisAnalysis->AddSigRegs( lowDMreg );
  ThisAnalysis->AddSigRegs( highDMreg );
  ThisAnalysis->AddSigRegs( inclusive );


  if( argument=="all" || argument=="output" || argument=="out" || argument=="table" || argument=="tables" ) makeTables( ThisAnalysis );
  if( argument=="all" || argument=="output" || argument=="out" || argument=="plots" || argument=="stacks" ) makeStack(  ThisAnalysis );

  return 0;
}
