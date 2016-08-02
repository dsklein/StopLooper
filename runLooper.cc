#include "analysis.h"
#include "sample.h"
#include "runLooper.h"
#include "sigRegion.h"
#include "sigRegion.cc"

#include "TChain.h"
#include "TFile.h"
#include "TString.h"

#include "CMS3.h"


// Global variables, for use in defining signal regions
extern bool j1_isBtag;
extern double j1pt;

bool j1_isBtag;
double j1pt;

// Help with program options
void printHelp() {
  std::cout << "\nUsage: ./runLooper [arg]\n" << std::endl;

  std::cout << "Takes zero or more arguments from the following list:" << std::endl;
  std::cout << "[blank]     equivalent to 'all'" << std::endl;
  std::cout << "help        show this message" << std::endl;
  std::cout << "all         run ScanChain, makeTables, makeStack, and makeDataCards" << std::endl;
  std::cout << "scan        run ScanChain only" << std::endl;
  std::cout << "plots       run makeStack only" << std::endl;
  std::cout << "tables      run makeTables only" << std::endl;
  std::cout << "cards       run makeDataCards only" << std::endl;
  std::cout << "estimate    run makeLostLepEstimate only" << std::endl;
  std::cout << "output      run makeStack, makeTables, and makeDataCards only" << std::endl;
}

///// MAIN PROGRAM /////
int main( int argc, char* argv[] ) {


  TString bkgPath = "/nfs-6/userdata/mliu/onelepbabies/V80_7p65_v2/";
  TString sigPath = "/hadoop/cms/store/user/haweber/condor/stop1l_2016/stop_babies_V080009_signal_norm_v2/merged_files/";
  TString dataPath = "/hadoop/cms/store/user/jgwood/condor/stop_1l_babies/stop_babies__CMS3_V080005__BabyMaker_V0800X_v1__20160612/merged_files/";


  ////////////////////////////////////////////////////////////////
  // Parse command-line argument(s)

  std::vector<TString> arguments;
  if( argc==1 )                       arguments.push_back( "all" );
  else  for( int i=1; i<argc; i++ )  arguments.push_back( TString(argv[i]) );

  bool needshelp   = false;
  bool runlooper   = false;
  bool runstacks   = false;
  bool runtables   = false;
  bool runcards    = false;
  bool runlostlep  = false;
  bool runestimate = false;


  for( TString arg : arguments ) {
	if(      arg=="help"  || arg=="h" ) needshelp = true;
	else if( arg=="scan"  || arg=="loop"  || arg=="scanchain" ) runlooper = true;
	else if( arg=="plot"  || arg=="plots" || arg=="stack" || arg=="stacks" ) runstacks = true;
	else if( arg=="table" || arg=="tables" ) runtables = true;
	else if( arg=="cards" || arg=="card"  || arg=="datacards" || arg=="datacard" ) runcards = true;
	else if( arg=="lostlep" || arg=="lost" || arg=="ll" ) runlostlep = true;
	else if( arg=="estimate" || arg=="est" || arg=="bkg" ) runestimate = true;
	else if( arg=="out"   || arg=="output" ) {
	  runstacks = true;
	  runtables = true;
	  runcards  = true;
	}
	else if( arg=="all") {
	  runlooper = true;
	  runstacks = true;
	  runtables = true;
	  runcards  = true;
	  runlostlep = true;
	  runestimate = true;
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
  // Make "analysis" objects out of "samples"

  analysis* srAnalysis = new analysis( 2.07, "plots.root" );
  analysis* crLostLep  = new analysis( 2.07, "plotsLL.root" );

  //                           new sample( "Label",  "Display name",    TColor,    sampleType )

  sample* data    =            new sample( "data",    "Data",           kBlack,    sample::kData );
  sample* signal  = srAnalysis->AddSample( "signal",  "T2tt",           kBlue+3,   sample::kSignal );

  sample* tt2l    = srAnalysis->AddSample( "tt2l", "$t\\bar{t} \\rightarrow 2l$", "t#bar{t} #rightarrow 2l", kCyan-3,   sample::kBackground );
  sample* tt1l    = srAnalysis->AddSample( "tt1l", "$t\\bar{t} \\rightarrow 1l$", "t#bar{t} #rightarrow 1l", kRed-7,    sample::kBackground );
  sample* singtop = srAnalysis->AddSample( "singletop", "Single Top",   kGreen-4,  sample::kBackground );
  sample* wjets   = srAnalysis->AddSample( "wjets",   "W+Jets",         kOrange-2, sample::kBackground );
  sample* dy      = srAnalysis->AddSample( "dy",      "Drell-Yan",      kRed+2,    sample::kBackground );
  sample* rare    = srAnalysis->AddSample( "rare",    "Rare",           kMagenta-5,sample::kBackground );

  // srAnalysis->AddSample( data );   // Uncomment this line to unblind
  crLostLep->AddSample( data );
  crLostLep->AddSample( signal );
  crLostLep->AddSample( tt2l );
  crLostLep->AddSample( tt1l );
  crLostLep->AddSample( singtop );
  crLostLep->AddSample( wjets );
  crLostLep->AddSample( dy );
  crLostLep->AddSample( rare );


  //////////////////////////////////////////////
  // Make selection and signal region objects
  selection<float> MET_250_350( (*tas::pfmet), 250., 350. );
  selection<float> MET_350_450( (*tas::pfmet), 350., 450. );
  selection<float> MET_450_550( (*tas::pfmet), 450., 550. );
  selection<float> MET_550_650( (*tas::pfmet), 550., 650. );
  selection<float> MET_350_inf( (*tas::pfmet), 350., 9999999. );
  selection<float> MET_450_inf( (*tas::pfmet), 450., 9999999. );
  selection<float> MET_550_inf( (*tas::pfmet), 550., 9999999. );
  selection<float> MET_650_inf( (*tas::pfmet), 650., 9999999. );

  selection<int> nJetsEq2( (*tas::ngoodjets), 2 );
  selection<int> nJetsEq3( (*tas::ngoodjets), 3 );
  selection<int> nJetsGe4( (*tas::ngoodjets), 4, 9999999 );
  selection<int> nJetsGe5( (*tas::ngoodjets), 5, 9999999 );

  selection<float>  lowMT2W( (*tas::MT2W),   0., 200.     );
  selection<float> highMT2W( (*tas::MT2W), 200., 9999999. );

  selection<float> modTop( (*tas::topnessMod), 6.4, 999999. );

  selection<double> j1Pt200( &j1pt, 200., 999999. );
  selection<bool>   j1NoTag( &j1_isBtag, false );


  sigRegion compr250(    "compr250",    "2 jets, modTop, MET 250-350" );
  sigRegion compr350(    "compr350",    "2 jets, modTop, MET 350-450" );
  sigRegion compr450(    "compr450",    "2 jets, modTop, MET 450+" );
  sigRegion boost250(    "boost250",    "3 jets, high MT2W, MET 250-350" );
  sigRegion boost350(    "boost350",    "3 jets, high MT2W, MET 350-450" );
  sigRegion boost450(    "boost450",    "3 jets, high MT2W, MET 450-550" );
  sigRegion boost550(    "boost550",    "3 jets, high MT2W, MET 550+" );
  sigRegion low250(      "low250",      "4+ jets, low MT2W, MET 250-350" );
  sigRegion low350(      "low350",      "4+ jets, low MT2W, MET 350-450" );
  sigRegion low450(      "low450",      "4+ jets, low MT2W, MET 450+" );
  sigRegion high250(     "high250",     "4+ jets, high MT2W, MET 250-350" );
  sigRegion high350(     "high350",     "4+ jets, high MT2W, MET 350-450" );
  sigRegion high450(     "high450",     "4+ jets, high MT2W, MET 450-550" );
  sigRegion high550(     "high550",     "4+ jets, high MT2W, MET 550-650" );
  sigRegion high650(     "high650",     "4+ jets, high MT2W, MET 650+" );
  sigRegion inclusive(   "inclusive",   "Inclusive" );
  sigRegion corridor250( "corridor250", "Corridor, low MET" );
  sigRegion corridor350( "corridor350", "Corridor, high MET" );

  compr250.AddSelections(    {&nJetsEq2, &MET_250_350, &modTop}   );
  compr350.AddSelections(    {&nJetsEq2, &MET_350_450, &modTop}   );
  compr450.AddSelections(    {&nJetsEq2, &MET_450_inf, &modTop}   );
  boost250.AddSelections(    {&nJetsEq3, &MET_250_350, &highMT2W} );
  boost350.AddSelections(    {&nJetsEq3, &MET_350_450, &highMT2W} );
  boost450.AddSelections(    {&nJetsEq3, &MET_450_550, &highMT2W} );
  boost550.AddSelections(    {&nJetsEq3, &MET_550_inf, &highMT2W} );
  low250.AddSelections(      {&nJetsGe4, &MET_250_350, &lowMT2W}  );
  low350.AddSelections(      {&nJetsGe4, &MET_350_450, &lowMT2W}  );
  low450.AddSelections(      {&nJetsGe4, &MET_450_inf, &lowMT2W}  );
  high250.AddSelections(     {&nJetsGe4, &MET_250_350, &highMT2W} );
  high350.AddSelections(     {&nJetsGe4, &MET_350_450, &highMT2W} );
  high450.AddSelections(     {&nJetsGe4, &MET_450_550, &highMT2W} );
  high550.AddSelections(     {&nJetsGe4, &MET_550_650, &highMT2W} );
  high650.AddSelections(     {&nJetsGe4, &MET_650_inf, &highMT2W} );
  corridor250.AddSelections( {&nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag} );
  corridor350.AddSelections( {&nJetsGe5, &MET_350_inf, &j1Pt200, &j1NoTag} );


  srAnalysis->AddSigRegs( {compr250, compr350, compr450} );
  srAnalysis->AddSigRegs( {boost250, boost350, boost450, boost550} );
  srAnalysis->AddSigRegs( {low250,  low350, low450} );
  srAnalysis->AddSigRegs( {high250, high350, high450, high550, high650} );
  srAnalysis->AddSigRegs( {inclusive} );
  srAnalysis->AddSigRegs( {corridor250, corridor350} );

  crLostLep->AddSigRegs( {compr250, compr350, compr450} );
  crLostLep->AddSigRegs( {boost250, boost350, boost450, boost550} );
  crLostLep->AddSigRegs( {low250,  low350, low450} );
  crLostLep->AddSigRegs( {high250, high350, high450, high550, high650} );
  crLostLep->AddSigRegs( {inclusive} );
  crLostLep->AddSigRegs( {corridor250, corridor350} );




  /////////////////////
  // Make chains

  if( runlooper || runlostlep ) {

	// Data samples
	data->AddFile( dataPath + "data_met_Run2016B_MINIAOD_PromptReco-v2.root" );
	data->AddFile( dataPath + "data_single_electron_Run2016B_MINIAOD_PromptReco-v2.root" );
	data->AddFile( dataPath + "data_single_muon_Run2016B_MINIAOD_PromptReco-v2.root" );

	// Signal sample(s)
	signal->AddFile( sigPath + "Signal_T2tt*.root" );

	// Background samples
	tt2l->AddFile( bkgPath + "ttbar_diLept_madgraph_pythia8_ext1*25ns*.root" );

	tt1l->AddFile( bkgPath + "ttbar_singleLeptFromT_madgraph_pythia8_*.root" );
	tt1l->AddFile( bkgPath + "ttbar_singleLeptFromTbar_madgraph_pythia8_ext1*.root" );

	wjets->AddFile( bkgPath + "WJetsToLNu_HT*.root" );

	dy->AddFile( bkgPath + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns.root" );
	dy->AddFile( bkgPath + "DYJetsToLL_m50_amcnlo_pythia8_25ns.root" );

	singtop->AddFile( bkgPath + "t_sch_4f_amcnlo_pythia8_25ns.root" );
	// singtop->AddFile( bkgPath + "t_tch_4f_powheg_pythia8_25ns.root" );
	// singtop->AddFile( bkgPath + "tbar_tch_4f_powheg_pythia8_25ns.root" );
	singtop->AddFile( bkgPath + "t_tW_5f_powheg_pythia8_25ns.root" );
	singtop->AddFile( bkgPath + "t_tbarW_5f_powheg_pythia8_25ns.root" );

	rare->AddFile( bkgPath + "TTWJetsToLNu_amcnlo_pythia8_25ns.root" );
	rare->AddFile( bkgPath + "TTWJetsToQQ_amcnlo_pythia8_25ns.root" );
	rare->AddFile( bkgPath + "TTZToLLNuNu_M-10_amcnlo_pythia8_25ns.root" );
	rare->AddFile( bkgPath + "TTZToQQ_amcnlo_pythia8_25ns.root" );
	// rare->AddFile( bkgPath + "tZq_ll_4f_amcnlo_pythia8_25ns.root" );
	// rare->AddFile( bkgPath + "tZq_nunu_4f_amcnlo_pythia8_25ns.root" );
	rare->AddFile( bkgPath + "WWTo2l2Nu_powheg_25ns.root" );
	rare->AddFile( bkgPath + "WWToLNuQQ_powheg_25ns.root" );
	rare->AddFile( bkgPath + "WZTo3LNu_powheg_pythia8_25ns.root" );
	rare->AddFile( bkgPath + "WZTo2L2Q_amcnlo_pythia8_25ns.root" );
	rare->AddFile( bkgPath + "WZTo1LNu2Q_amcnlo_pythia8_25ns.root" );
	// rare->AddFile( bkgPath + "ZZTo4L_powheg_pythia8_25ns.root" );
	rare->AddFile( bkgPath + "ZZTo2L2Q_amcnlo_pythia8_25ns.root" );
	rare->AddFile( bkgPath + "ZZTo2L2Nu_powheg_pythia8_25ns.root" );
	rare->AddFile( bkgPath + "ZZTo2Q2Nu_amcnlo_pythia8_25ns.root" );
  }

  ///////////////////////////////////////
  // Run ScanChain (signal regions)

  if( runlooper ) {

	// Reset output file
	TFile* outfile = new TFile( srAnalysis->GetFileName(), "RECREATE");
	outfile->Close();

	// Run ScanChain on all samples
	for( sample* mySample : srAnalysis->GetAllSamples() ) ScanChain( srAnalysis, mySample );

  }

  //////////////////////////////////////////////
  // Run lost lepton background estimate

  if( runlostlep ) {

	// Reset output file
	TFile* outfile = new TFile( crLostLep->GetFileName(), "RECREATE");
	outfile->Close();

	// Run lost lepton CR looper on all samples
	for( sample* mySample : crLostLep->GetAllSamples() ) looperCR2lep( crLostLep, mySample );

  }


  //////////////////////////////////
  // Make all the various outputs 

  if( runtables   ) makeTables( srAnalysis );
  if( runlostlep || runtables ) makeTables( crLostLep  );
  if( runstacks   ) makeStack( srAnalysis );
  if( runestimate ) makeLostLepEstimate( srAnalysis, crLostLep );
  if( runcards    ) makeDataCards( srAnalysis );

  return 0;
}
