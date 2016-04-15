#include "analysis.h"
#include "sample.h"
#include "runLooper.h"

#include "TChain.h"
#include "TFile.h"
#include "TString.h"

// Help with program options
void printHelp() {
  std::cout << "\nUsage: ./runLooper [arg]\n" << std::endl;

  std::cout << "Takes zero or more arguments from the following list:" << std::endl;
  std::cout << "[blank]     equivalent to 'all'" << std::endl;
  std::cout << "help        show this message" << std::endl;
  std::cout << "all         run ScanChain, makeTables, and makeStack" << std::endl;
  std::cout << "scan        run ScanChain only" << std::endl;
  std::cout << "plots       run makeStack only" << std::endl;
  std::cout << "tables      run makeTables only" << std::endl;
  std::cout << "output      run makeStack and makeTables only" << std::endl;
}

///// MAIN PROGRAM /////
int main( int argc, char* argv[] ) {


  TString bkgPath = "/hadoop/cms/store/user/jgwood/condor/stop_1l_babies/stop_babies__CMS3_V070411__BabyMaker_V0704X_v9__20160127/merged_files/Skims_SR__20160202/";

  TString sigPath = "/nfs-7/userdata/stopRun2/signalbabies/";


  ////////////////////////////////////////////////////////////////
  // Parse command-line argument(s)

  std::vector<TString> arguments;
  if( argc==1 )                       arguments.push_back( "all" );
  else  for( int i=1; i<argc; i++ )  arguments.push_back( TString(argv[i]) );

  bool needshelp = false;
  bool runlooper = false;
  bool runstacks = false;
  bool runtables = false;


  for( TString arg : arguments ) {
	if(      arg=="help"  || arg=="h" ) needshelp = true;
	else if( arg=="scan"  || arg=="loop"  || arg=="scanchain" ) runlooper = true;
	else if( arg=="plot"  || arg=="plots" || arg=="stack" || arg=="stacks" ) runstacks = true;
	else if( arg=="table" || arg=="tables" ) runtables = true;
	else if( arg=="out"   || arg=="output" ) {
	  runstacks = true;
	  runtables = true;
	}
	else if( arg=="all") {
	  runlooper = true;
	  runstacks = true;
	  runtables = true;
	}
	else {
	  std::cout << "Unrecognized option: " << arg << std::endl;
	  needshelp = true;
	}
  }

  // If the user inputs a wrong option, print the help and exit
  if( needshelp ) {
	printHelp();
	return 1;
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // Make "analysis" object out of "samples"

  analysis* ThisAnalysis = new analysis( 2.26 ); // Luminosity

  //                             new sample( "Label",  "Display name",    TColor,    sampleType )

  sample* stop700 = ThisAnalysis->AddSample( "stop700", "T2tt (700,50)",  kBlue+3,   sample::kSignal );
  sample* stop600 = ThisAnalysis->AddSample( "stop600", "T2tt (600,250)", kGreen+3,  sample::kSignal );
  sample* stop300 = ThisAnalysis->AddSample( "stop300", "T2tt (300,200)", kMagenta+3,sample::kSignal );
  sample* stop275 = ThisAnalysis->AddSample( "stop275", "T2tt (275,100)", kOrange+7, sample::kSignal );

  sample* tt2l    = ThisAnalysis->AddSample( "tt2l", "$t\\bar{t} \\rightarrow 2l$", "t#bar{t} #rightarrow 2l", kCyan-3,   sample::kBackground );
  sample* tt1l    = ThisAnalysis->AddSample( "tt1l", "$t\\bar{t} \\rightarrow 1l$", "t#bar{t} #rightarrow 1l", kRed-7,    sample::kBackground );
  sample* singtop = ThisAnalysis->AddSample( "singletop", "Single Top",   kGreen-4,  sample::kBackground );
  sample* wjets   = ThisAnalysis->AddSample( "wjets",   "W+Jets",         kOrange-2, sample::kBackground );
  sample* dy      = ThisAnalysis->AddSample( "dy",      "Drell-Yan",      kRed+2,    sample::kBackground );
  sample* rare    = ThisAnalysis->AddSample( "rare",    "Rare",           kMagenta-5,sample::kBackground );

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



  //////////////////////////////////////////////
  // Make chains and run ScanChain

  if( runlooper ) {

	// Signal samples

	stop700->AddFile( sigPath + "Signal_T2tt_700_50.root");
	stop600->AddFile( sigPath + "Signal_T2tt_600_250.root");
	stop300->AddFile( sigPath + "Signal_T2tt_300_200.root");
	stop275->AddFile( sigPath + "Signal_T2tt_275_100.root");


	// Background samples

	tt2l->AddFile( bkgPath + "ttbar_powheg_pythia8_25ns_skimmed.root" );

	tt1l->AddFile( bkgPath + "ttbar_powheg_pythia8_25ns_skimmed.root" );

	wjets->AddFile( bkgPath + "WJetsToLNu_HT100To200_madgraph_pythia8_25ns_skimmed.root" );
	wjets->AddFile( bkgPath + "WJetsToLNu_HT200To400_madgraph_pythia8_25ns_skimmed.root" );
	wjets->AddFile( bkgPath + "WJetsToLNu_HT400To600_madgraph_pythia8_25ns_skimmed.root" );
	wjets->AddFile( bkgPath + "WJetsToLNu_HT600ToInf_madgraph_pythia8_25ns_skimmed.root" );

	dy->AddFile( bkgPath + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns_skimmed.root" );
	dy->AddFile( bkgPath + "DYJetsToLL_m50_amcnlo_pythia8_25ns_skimmed.root" );

	singtop->AddFile( bkgPath + "t_sch_4f_amcnlo_pythia8_25ns_skimmed.root" );
	singtop->AddFile( bkgPath + "t_tch_4f_powheg_pythia8_25ns_skimmed.root" );
	singtop->AddFile( bkgPath + "tbar_tch_4f_powheg_pythia8_25ns_skimmed.root" );
	singtop->AddFile( bkgPath + "t_tW_5f_powheg_pythia8_25ns_skimmed.root" );
	singtop->AddFile( bkgPath + "t_tbarW_5f_powheg_pythia8_25ns_skimmed.root" );

	rare->AddFile( bkgPath + "TTWJetsToLNu_amcnlo_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "TTWJetsToQQ_amcnlo_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "TTZToLLNuNu_M-10_amcnlo_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "TTZToQQ_amcnlo_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "tZq_ll_4f_amcnlo_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "tZq_nunu_4f_amcnlo_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "WWTo2l2Nu_powheg_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "WWToLNuQQ_powheg_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "WZTo3LNu_powheg_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "WZTo2L2Q_amcnlo_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "WZTo1LNu2Q_amcnlo_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "ZZTo4L_powheg_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "ZZTo2L2Q_amcnlo_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "ZZTo2L2Nu_powheg_pythia8_25ns_skimmed.root" );
	rare->AddFile( bkgPath + "ZZTo2Q2Nu_amcnlo_pythia8_25ns_skimmed.root" );


	// Reset output file and run ScanChain on all samples
	TFile* outfile = new TFile("plots.root", "RECREATE");
	outfile->Close();

	ScanChain( ThisAnalysis, stop700 );
	ScanChain( ThisAnalysis, stop600 );
	ScanChain( ThisAnalysis, stop300 );
	ScanChain( ThisAnalysis, stop275 );
	ScanChain( ThisAnalysis, tt2l    );
	ScanChain( ThisAnalysis, tt1l    );
	ScanChain( ThisAnalysis, singtop );
	ScanChain( ThisAnalysis, wjets   );
	ScanChain( ThisAnalysis, dy      );
	ScanChain( ThisAnalysis, rare    );
  }


  /////////////////////////////////////////////////////
  // Make stacked histograms and/or yield tables

  if( runtables ) makeTables( ThisAnalysis );
  if( runstacks ) makeStack(  ThisAnalysis );

  return 0;
}
