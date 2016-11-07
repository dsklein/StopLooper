#include "analysis.h"
#include "sample.h"
#include "runLooper.h"
#include "sigRegion.h"
#include "sfHelper.h"
#include "systematic.h"

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
	std::cout << "all         run all of the components listed below" << std::endl;
	std::cout << "scan        run ScanChain only" << std::endl;
	std::cout << "lostlep     run lost lepton looper only" << std::endl;
	std::cout << "1lw         run 1l-from-W background looper only" << std::endl;
	std::cout << "jes         run loopers to calculate JES systematic variations" << std::endl;
	std::cout << "plots       run makeStack only" << std::endl;
	std::cout << "tables      run makeTables only" << std::endl;
	std::cout << "estimate    run data-driven background estimates only" << std::endl;
	std::cout << "cards       run makeDataCards only" << std::endl;
	std::cout << "output      run makeStack, makeTables, and makeDataCards only" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
//------------------------------ MAIN PROGRAM --------------------------------//

int main( int argc, char* argv[] ) {


	// Parse command-line argument(s)

	std::vector<TString> arguments;
	if( argc==1 )                      arguments.push_back( "all" );
	else  for( int i=1; i<argc; i++ )  arguments.push_back( TString(argv[i]) );

	bool needshelp   = false;
	bool runlooper   = false;
	bool runlostlep  = false;
	bool run1lw      = false;
	bool runjes      = false;
	bool runstacks   = false;
	bool runtables   = false;
	bool runcards    = false;
	bool runestimate = false;


	for( TString arg : arguments ) {
		if(      arg=="help"  || arg=="h" ) needshelp = true;
		else if( arg=="scan"  || arg=="loop"  || arg=="scanchain" ) runlooper = true;
		else if( arg=="lostlep" || arg=="lost" || arg=="ll" ) runlostlep = true;
		else if( arg=="1lw" || arg=="wjets" || arg=="0bjets" ) run1lw = true;
		else if( arg=="jes" || arg=="JES" ) runjes = true;
		else if( arg=="plot"  || arg=="plots" || arg=="stack" || arg=="stacks" ) runstacks = true;
		else if( arg=="table" || arg=="tables" ) runtables = true;
		else if( arg=="cards" || arg=="card"  || arg=="datacards" || arg=="datacard" ) runcards = true;
		else if( arg=="estimate" || arg=="est" || arg=="bkg" ) runestimate = true;
		else if( arg=="out"   || arg=="output" ) {
			runstacks = true;
			runtables = true;
			runcards  = true;
		}
		else if( arg=="all") {
			runlooper = true;
			runlostlep = true;
			run1lw  = true;
			runjes    = true;
			runstacks = true;
			runtables = true;
			runcards  = true;
			runestimate = true;
		}
		else {
			std::cout << "Unrecognized option: " << arg << std::endl;
			needshelp = true;
		}
	}

	// If the user inputs a wrong option or asks for help, print the list of allowed options and exit
	if( needshelp ) {
		printHelp();
		return 1;
	}


	////////////////////////////////////////////////////////////////////////////////////////
	// Define the "sample" and "analysis" objects that will do much of our bookkeeping

	//                     new analysis( lumi, "histogram storage file", "systematic storage file" )
	analysis* srAnalysis = new analysis( 12.9, "plots.root", "systVariations.root" );
	analysis* crLostLep  = new analysis( 12.9, "plotsLL.root", "systVariationsLL.root" );
	analysis* cr0bjets;     // Will be defined later
	analysis* sr_jesup   = new analysis( 12.9, "jes_sr.root", "jes_sr.root" );
	analysis* sr_jesdn   = new analysis( 12.9, "jes_sr.root", "jes_sr.root" );
	analysis* cr2l_jesup = new analysis( 12.9, "jes_cr2l.root", "jes_cr2l.root" );
	analysis* cr2l_jesdn = new analysis( 12.9, "jes_cr2l.root", "jes_cr2l.root" );
	analysis* cr0b_jesup;
	analysis* cr0b_jesdn;

	//                new sample( "Label",  "Display name(s)", TColor,    sampleType )
	sample* data    = new sample( "data",    "Data",           kBlack,    sample::kData );
	sample* signal  = new sample( "signal",  "T2tt",           kBlue+3,   sample::kSignal );
	sample* tt2l    = new sample( "tt2l", "$t\\bar{t} \\rightarrow 2l$", "t#bar{t} #rightarrow 2l", kCyan-3,   sample::kBackground );
	sample* tt1l    = new sample( "tt1l", "$t\\bar{t} \\rightarrow 1l$", "t#bar{t} #rightarrow 1l", kRed-7,    sample::kBackground );
	sample* singtop = new sample( "singletop", "Single Top",   kGreen-4,  sample::kBackground );
	sample* wjets   = new sample( "wjets",   "W+Jets",         kOrange-2, sample::kBackground );
	sample* dy      = new sample( "dy",      "Drell-Yan",      kRed+2,    sample::kBackground );
	sample* rare    = new sample( "rare",    "Rare",           kMagenta-5,sample::kBackground );

	sample* signal_jesup  = new sample( "signal_jesup",  "T2tt",           kBlue+3,   sample::kSignal );
	sample* tt2l_jesup    = new sample( "tt2l_jesup", "$t\\bar{t} \\rightarrow 2l$", "t#bar{t} #rightarrow 2l", kCyan-3,   sample::kBackground );
	sample* tt1l_jesup    = new sample( "tt1l_jesup", "$t\\bar{t} \\rightarrow 1l$", "t#bar{t} #rightarrow 1l", kRed-7,    sample::kBackground );
	sample* singtop_jesup = new sample( "singletop_jesup", "Single Top",   kGreen-4,  sample::kBackground );
	sample* wjets_jesup   = new sample( "wjets_jesup",   "W+Jets",         kOrange-2, sample::kBackground );
	sample* dy_jesup      = new sample( "dy_jesup",      "Drell-Yan",      kRed+2,    sample::kBackground );
	sample* rare_jesup    = new sample( "rare_jesup",    "Rare",           kMagenta-5,sample::kBackground );

	sample* signal_jesdn  = new sample( "signal_jesdn",  "T2tt",           kBlue+3,   sample::kSignal );
	sample* tt2l_jesdn    = new sample( "tt2l_jesdn", "$t\\bar{t} \\rightarrow 2l$", "t#bar{t} #rightarrow 2l", kCyan-3,   sample::kBackground );
	sample* tt1l_jesdn    = new sample( "tt1l_jesdn", "$t\\bar{t} \\rightarrow 1l$", "t#bar{t} #rightarrow 1l", kRed-7,    sample::kBackground );
	sample* singtop_jesdn = new sample( "singletop_jesdn", "Single Top",   kGreen-4,  sample::kBackground );
	sample* wjets_jesdn   = new sample( "wjets_jesdn",   "W+Jets",         kOrange-2, sample::kBackground );
	sample* dy_jesdn      = new sample( "dy_jesdn",      "Drell-Yan",      kRed+2,    sample::kBackground );
	sample* rare_jesdn    = new sample( "rare_jesdn",    "Rare",           kMagenta-5,sample::kBackground );

	// srAnalysis->AddSample( data );   // Uncomment this line to unblind
	crLostLep->AddSample( data );
	srAnalysis->AddSample( signal );   crLostLep->AddSample( signal );
	srAnalysis->AddSample( tt2l );     crLostLep->AddSample( tt2l );
	srAnalysis->AddSample( tt1l );     crLostLep->AddSample( tt1l );
	srAnalysis->AddSample( singtop );  crLostLep->AddSample( singtop );
	srAnalysis->AddSample( wjets );    crLostLep->AddSample( wjets );
	srAnalysis->AddSample( dy );       crLostLep->AddSample( dy );
	srAnalysis->AddSample( rare );     crLostLep->AddSample( rare );

	sr_jesup->AddSample( signal_jesup );   cr2l_jesup->AddSample( signal_jesup );
	sr_jesup->AddSample( tt2l_jesup );     cr2l_jesup->AddSample( tt2l_jesup );
	sr_jesup->AddSample( tt1l_jesup );     cr2l_jesup->AddSample( tt1l_jesup );
	sr_jesup->AddSample( singtop_jesup );  cr2l_jesup->AddSample( singtop_jesup );
	sr_jesup->AddSample( wjets_jesup );    cr2l_jesup->AddSample( wjets_jesup );
	sr_jesup->AddSample( dy_jesup );       cr2l_jesup->AddSample( dy_jesup );
	sr_jesup->AddSample( rare_jesup );     cr2l_jesup->AddSample( rare_jesup );

	sr_jesdn->AddSample( signal_jesdn );   cr2l_jesdn->AddSample( signal_jesdn );
	sr_jesdn->AddSample( tt2l_jesdn );     cr2l_jesdn->AddSample( tt2l_jesdn );
	sr_jesdn->AddSample( tt1l_jesdn );     cr2l_jesdn->AddSample( tt1l_jesdn );
	sr_jesdn->AddSample( singtop_jesdn );  cr2l_jesdn->AddSample( singtop_jesdn );
	sr_jesdn->AddSample( wjets_jesdn );    cr2l_jesdn->AddSample( wjets_jesdn );
	sr_jesdn->AddSample( dy_jesdn );       cr2l_jesdn->AddSample( dy_jesdn );
	sr_jesdn->AddSample( rare_jesdn );     cr2l_jesdn->AddSample( rare_jesdn );

	/////////////////////////////////////////////////////////////////////////
	// Create "systematic" objects to store all our systematic variations

	//                         "Label",    variation direction, function that provides the variation
	systematic jesup(          "JES",      systematic::kSkipUp,   NULL );
	systematic jesdn(          "JES",      systematic::kSkipDown, NULL );
	systematic lepSFup(        "lepSF",    systematic::kUp,    (*sfhelp::LepSFUp) );
	systematic lepSFdn(        "lepSF",    systematic::kDown,  (*sfhelp::LepSFDown) );
	systematic btagHFup(       "btagHF",   systematic::kUp,    (*sfhelp::BtagHeavyUp) );
	systematic btagHFdn(       "btagHF",   systematic::kDown,  (*sfhelp::BtagHeavyDown) );
	systematic btagLFup(       "btagLF",   systematic::kUp,    (*sfhelp::BtagLightUp) );
	systematic btagLFdn(       "btagLF",   systematic::kDown,  (*sfhelp::BtagLightDown) );
	systematic qSquaredup(     "qSquared", systematic::kUp,    (*sfhelp::QSquaredUp) );
	systematic qSquareddn(     "qSquared", systematic::kDown,  (*sfhelp::QSquaredDown) );
	systematic alphaSup(       "alphaS",   systematic::kUp,    (*sfhelp::AlphaSUp) );
	systematic alphaSdn(       "alphaS",   systematic::kDown,  (*sfhelp::AlphaSDown) );
	systematic eff2lup(        "cr2ltrig", systematic::kUp,    (*sfhelp::Trig2lUp) );
	systematic eff2ldn(        "cr2ltrig", systematic::kDown,  (*sfhelp::Trig2lDown) );
	systematic eff2lup_dummy(  "cr2ltrig", systematic::kUp,    (*sfhelp::Unity) );
	systematic eff2ldn_dummy(  "cr2ltrig", systematic::kDown,  (*sfhelp::Unity) );
	systematic metresup(       "METres",   systematic::kUp,    (*sfhelp::MetResUp) );
	systematic metresdn(       "METres",   systematic::kDown,  (*sfhelp::MetResDown) );
	systematic topptup(        "topSysPt", systematic::kUp,    (*sfhelp::TopSystPtUp) );
	systematic topptdn(        "topSysPt", systematic::kDown,  (*sfhelp::TopSystPtDown) );
	systematic contam1lwup(    "contam",   systematic::kUp,    (*sfhelp::Contam1lwUp) );
	systematic contam1lwdn(    "contam",   systematic::kDown,  (*sfhelp::Contam1lwDown) );
	systematic contamup_dummy( "contam",   systematic::kUp,    (*sfhelp::Unity) );
	systematic contamdn_dummy( "contam",   systematic::kDown,  (*sfhelp::Unity) );

	srAnalysis->AddSystematics( {&jesup, &jesdn, &lepSFup, &lepSFdn, &btagHFup, &btagHFdn, &btagLFup, &btagLFdn, &qSquaredup, &qSquareddn, &alphaSup, &alphaSdn} );
	srAnalysis->AddSystematics( {&eff2lup_dummy, &eff2ldn_dummy, &metresup, &metresdn, &topptup, &topptdn, &contamup_dummy, &contamdn_dummy } );
	crLostLep->AddSystematics(  {&jesup, &jesdn, &lepSFup, &lepSFdn, &btagHFup, &btagHFdn, &btagLFup, &btagLFdn, &qSquaredup, &qSquareddn} );
	crLostLep->AddSystematics(  {&alphaSup, &alphaSdn, &eff2lup, &eff2ldn, &metresup, &metresdn, &topptup, &topptdn } );

	// A sneaky trick to make JES systematics work with existing code
	systematic jesup_dummy( "JES", systematic::kUp,   (*sfhelp::Unity) );
	systematic jesdn_dummy( "JES", systematic::kDown, (*sfhelp::Unity) );
	sr_jesup->AddSystematics( {&jesup_dummy} );
	sr_jesdn->AddSystematics( {&jesdn_dummy} );
	cr2l_jesup->AddSystematics( {&jesup_dummy} );
	cr2l_jesdn->AddSystematics( {&jesdn_dummy} );


	/////////////////////////////////////////////////////////////////////////////////////
	// Create the objects that will define our signal and control regions

	// Create "selection"s - objects that encode a cut on a baby branch or global variable
	// selection<type>  obj_name( (cutVariable), [minvalue, maxvalue] OR [equal value] )
	selection<float> MET_250_350( (*tas::pfmet), 250., 350. );
	selection<float> MET_350_450( (*tas::pfmet), 350., 450. );
	selection<float> MET_450_550( (*tas::pfmet), 450., 550. );
	selection<float> MET_550_650( (*tas::pfmet), 550., 650. );
	selection<float> MET_450_inf( (*tas::pfmet), 450., 9999999. ); // MET bins for the signal regions
	selection<float> MET_550_inf( (*tas::pfmet), 550., 9999999. );
	selection<float> MET_650_inf( (*tas::pfmet), 650., 9999999. );

	selection<float> CR_MET_250_350( (*tas::pfmet_rl), 250., 350. );
	selection<float> CR_MET_350_450( (*tas::pfmet_rl), 350., 450. );
	selection<float> CR_MET_450_550( (*tas::pfmet_rl), 450., 550. );
	selection<float> CR_MET_550_650( (*tas::pfmet_rl), 550., 650. );
	selection<float> CR_MET_450_inf( (*tas::pfmet_rl), 450., 9999999. ); // (MET+lep2) bins for the 2-lep control regions
	selection<float> CR_MET_550_inf( (*tas::pfmet_rl), 550., 9999999. );
	selection<float> CR_MET_650_inf( (*tas::pfmet_rl), 650., 9999999. );

	selection<int> nJetsEq2( (*tas::ngoodjets), 2 );
	selection<int> nJetsEq3( (*tas::ngoodjets), 3 );
	selection<int> nJetsGe4( (*tas::ngoodjets), 4, 9999999 ); // NJets bins
	selection<int> nJetsGe5( (*tas::ngoodjets), 5, 9999999 );

	selection<float>  lowMT2W(    (*tas::MT2W),      0., 200.     );
	selection<float> highMT2W(    (*tas::MT2W),    200., 9999999. );
	selection<float>  CR_lowMT2W( (*tas::MT2W_rl),   0., 200.     ); // MT2W bins
	selection<float> CR_highMT2W( (*tas::MT2W_rl), 200., 9999999. );

	selection<float> modTop(    (*tas::topnessMod),    6.4, 999999. ); // Modified topness for compressed T2tb regions
	selection<float> CR_modTop( (*tas::topnessMod_rl), 6.4, 999999. );

	selection<double> j1Pt200( &j1pt, 200., 999999. ); // Special selections for the corridor regions
	selection<bool>   j1NoTag( &j1_isBtag, false );


	// Create the "sigRegion" objects that will store the definitions of our signal/control regions
	// sigRegion objName(  "label",       "Nice name(s) for plots/tables",  {selections that define the region} )
	sigRegion compr250(    "compr250",    "2 jets, modTop, MET 250-350",     {&nJetsEq2, &MET_250_350, &modTop}   );
	sigRegion compr350(    "compr350",    "2 jets, modTop, MET 350-450",     {&nJetsEq2, &MET_350_450, &modTop}   );
	sigRegion compr450(    "compr450",    "2 jets, modTop, MET 450+",        {&nJetsEq2, &MET_450_inf, &modTop}   );
	sigRegion boost250(    "boost250",    "3 jets, high MT2W, MET 250-350",  {&nJetsEq3, &MET_250_350, &highMT2W} );
	sigRegion boost350(    "boost350",    "3 jets, high MT2W, MET 350-450",  {&nJetsEq3, &MET_350_450, &highMT2W} );
	sigRegion boost450(    "boost450",    "3 jets, high MT2W, MET 450-550",  {&nJetsEq3, &MET_450_550, &highMT2W} );
	sigRegion boost550(    "boost550",    "3 jets, high MT2W, MET 550+",     {&nJetsEq3, &MET_550_inf, &highMT2W} );
	sigRegion low250(      "low250",      "4+ jets, low MT2W, MET 250-350",  {&nJetsGe4, &MET_250_350, &lowMT2W}  );
	sigRegion low350(      "low350",      "4+ jets, low MT2W, MET 350-450",  {&nJetsGe4, &MET_350_450, &lowMT2W}  );
	sigRegion low450(      "low450",      "4+ jets, low MT2W, MET 450+",     {&nJetsGe4, &MET_450_inf, &lowMT2W}  );
	sigRegion high250(     "high250",     "4+ jets, high MT2W, MET 250-350", {&nJetsGe4, &MET_250_350, &highMT2W} );
	sigRegion high350(     "high350",     "4+ jets, high MT2W, MET 350-450", {&nJetsGe4, &MET_350_450, &highMT2W} );
	sigRegion high450(     "high450",     "4+ jets, high MT2W, MET 450-550", {&nJetsGe4, &MET_450_550, &highMT2W} );
	sigRegion high550(     "high550",     "4+ jets, high MT2W, MET 550-650", {&nJetsGe4, &MET_550_650, &highMT2W} );
	sigRegion high650(     "high650",     "4+ jets, high MT2W, MET 650+",    {&nJetsGe4, &MET_650_inf, &highMT2W} );
	sigRegion inclusive(   "inclusive",   "Inclusive" );
	sigRegion corridor250( "corridor250", "Corridor, low MET",  {&nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag} );
	sigRegion corridor350( "corridor350", "Corridor, mid MET",  {&nJetsGe5, &MET_350_450, &j1Pt200, &j1NoTag} );
	sigRegion corridor450( "corridor450", "Corridor, high MET", {&nJetsGe5, &MET_450_inf, &j1Pt200, &j1NoTag} );

	sigRegion compr250CR(    "compr250CR",    "CR 2 jets, modTop, MET 250-350",     {&nJetsEq2, &CR_MET_250_350, &CR_modTop}   );
	sigRegion compr350CR(    "compr350CR",    "CR 2 jets, modTop, MET 350-450",     {&nJetsEq2, &CR_MET_350_450, &CR_modTop}   );
	sigRegion compr450CR(    "compr450CR",    "CR 2 jets, modTop, MET 450+",        {&nJetsEq2, &CR_MET_450_inf, &CR_modTop}   );
	sigRegion boost250CR(    "boost250CR",    "CR 3 jets, high MT2W, MET 250-350",  {&nJetsEq3, &CR_MET_250_350, &CR_highMT2W} );
	sigRegion boost350CR(    "boost350CR",    "CR 3 jets, high MT2W, MET 350-450",  {&nJetsEq3, &CR_MET_350_450, &CR_highMT2W} );
	sigRegion boost450CR(    "boost450CR",    "CR 3 jets, high MT2W, MET 450-550",  {&nJetsEq3, &CR_MET_450_550, &CR_highMT2W} );
	sigRegion boost550CR(    "boost550CR",    "CR 3 jets, high MT2W, MET 550+",     {&nJetsEq3, &CR_MET_550_inf, &CR_highMT2W} );
	sigRegion low250CR(      "low250CR",      "CR 4+ jets, low MT2W, MET 250-350",  {&nJetsGe4, &CR_MET_250_350, &CR_lowMT2W}  );
	sigRegion low350CR(      "low350CR",      "CR 4+ jets, low MT2W, MET 350-450",  {&nJetsGe4, &CR_MET_350_450, &CR_lowMT2W}  );
	sigRegion low450CR(      "low450CR",      "CR 4+ jets, low MT2W, MET 450+",     {&nJetsGe4, &CR_MET_450_inf, &CR_lowMT2W}  );
	sigRegion high250CR(     "high250CR",     "CR 4+ jets, high MT2W, MET 250-350", {&nJetsGe4, &CR_MET_250_350, &CR_highMT2W} );
	sigRegion high350CR(     "high350CR",     "CR 4+ jets, high MT2W, MET 350-450", {&nJetsGe4, &CR_MET_350_450, &CR_highMT2W} );
	sigRegion high450CR(     "high450CR",     "CR 4+ jets, high MT2W, MET 450-550", {&nJetsGe4, &CR_MET_450_550, &CR_highMT2W} );
	sigRegion high550CR(     "high550CR",     "CR 4+ jets, high MT2W, MET 550-650", {&nJetsGe4, &CR_MET_550_650, &CR_highMT2W} );
	sigRegion high650CR(     "high650CR",     "CR 4+ jets, high MT2W, MET 650+",    {&nJetsGe4, &CR_MET_650_inf, &CR_highMT2W} );
	// sigRegion inclusive(   "inclusive",   "Inclusive" );
	sigRegion corridor250CR( "corridor250CR", "CR Corridor, low MET",  {&nJetsGe5, &CR_MET_250_350, &j1Pt200, &j1NoTag} );
	sigRegion corridor350CR( "corridor350CR", "CR Corridor, mid MET",  {&nJetsGe5, &CR_MET_350_450, &j1Pt200, &j1NoTag} );
	sigRegion corridor450CR( "corridor450CR", "CR Corridor, high MET", {&nJetsGe5, &CR_MET_450_inf, &j1Pt200, &j1NoTag} );


	// Finally, store all these signal/control regions in our "analysis" objects.
	// Each {vector of "sigRegions"} will give rise to its own yield table, so structure matters here!
	srAnalysis->AddSigRegs( {&compr250, &compr350, &compr450} );
	srAnalysis->AddSigRegs( {&boost250, &boost350, &boost450, &boost550} );
	srAnalysis->AddSigRegs( {&low250, & low350, &low450} );
	srAnalysis->AddSigRegs( {&high250, &high350, &high450, &high550, &high650} );
	srAnalysis->AddSigRegs( {&inclusive} );
	srAnalysis->AddSigRegs( {&corridor250, &corridor350, &corridor450} );

	crLostLep->AddSigRegs( {&compr250CR, &compr350CR, &compr450CR} );
	crLostLep->AddSigRegs( {&boost250CR, &boost350CR, &boost450CR, &boost550CR} );
	crLostLep->AddSigRegs( {&low250CR, & low350CR, &low450CR} );
	crLostLep->AddSigRegs( {&high250CR, &high350CR, &high450CR, &high550CR, &high650CR} );
	crLostLep->AddSigRegs( {&inclusive} );
	crLostLep->AddSigRegs( {&corridor250CR, &corridor350CR, &corridor450CR} );

	sr_jesup->AddSigRegs( {&compr250, &compr350, &compr450} );
	sr_jesup->AddSigRegs( {&boost250, &boost350, &boost450, &boost550} );
	sr_jesup->AddSigRegs( {&low250, & low350, &low450} );
	sr_jesup->AddSigRegs( {&high250, &high350, &high450, &high550, &high650} );
	sr_jesup->AddSigRegs( {&inclusive} );
	sr_jesup->AddSigRegs( {&corridor250, &corridor350, &corridor450} );

	cr2l_jesup->AddSigRegs( {&compr250CR, &compr350CR, &compr450CR} );
	cr2l_jesup->AddSigRegs( {&boost250CR, &boost350CR, &boost450CR, &boost550CR} );
	cr2l_jesup->AddSigRegs( {&low250CR, & low350CR, &low450CR} );
	cr2l_jesup->AddSigRegs( {&high250CR, &high350CR, &high450CR, &high550CR, &high650CR} );
	cr2l_jesup->AddSigRegs( {&inclusive} );
	cr2l_jesup->AddSigRegs( {&corridor250CR, &corridor350CR, &corridor450CR} );

	sr_jesdn->AddSigRegs( {&compr250, &compr350, &compr450} );
	sr_jesdn->AddSigRegs( {&boost250, &boost350, &boost450, &boost550} );
	sr_jesdn->AddSigRegs( {&low250, & low350, &low450} );
	sr_jesdn->AddSigRegs( {&high250, &high350, &high450, &high550, &high650} );
	sr_jesdn->AddSigRegs( {&inclusive} );
	sr_jesdn->AddSigRegs( {&corridor250, &corridor350, &corridor450} );

	cr2l_jesdn->AddSigRegs( {&compr250CR, &compr350CR, &compr450CR} );
	cr2l_jesdn->AddSigRegs( {&boost250CR, &boost350CR, &boost450CR, &boost550CR} );
	cr2l_jesdn->AddSigRegs( {&low250CR, & low350CR, &low450CR} );
	cr2l_jesdn->AddSigRegs( {&high250CR, &high350CR, &high450CR, &high550CR, &high650CR} );
	cr2l_jesdn->AddSigRegs( {&inclusive} );
	cr2l_jesdn->AddSigRegs( {&corridor250CR, &corridor350CR, &corridor450CR} );



	//////////////////////////////////////////////////////////////////////////////////////////////
	// For each "sample" object defined earlier, chain up the baby files that make up that sample

	TString sigPath = "/hadoop/cms/store/user/haweber/condor/stop1l_2016/stop_babies_V080009_signal_norm_v2/merged_files/";
	TString bkgPath = "/nfs-7/userdata/stopRun2/analysis2016__SUS-16-028__12p9fb/stopBabies__v8.0.x_v8__20160729/Nominal/";
	TString dataPath = "/nfs-7/userdata/stopRun2/analysis2016__SUS-16-028__12p9fb/stopBabies__v8.0.x_v8__20160729/Nominal/";

	TString sigPath_jesup = "/hadoop/cms/store/user/haweber/condor/stop1l_2016/stop_babies_V080009_signal_JESup_v2/merged_files/";
	TString bkgPath_jesup = "/nfs-7/userdata/stopRun2/analysis2016__SUS-16-028__12p9fb/stopBabies__v8.0.x_v8__20160729/JESup/";
	TString sigPath_jesdn = "/hadoop/cms/store/user/haweber/condor/stop1l_2016/stop_babies_V080009_signal_JESdown_v2/merged_files/";
	TString bkgPath_jesdn = "/nfs-7/userdata/stopRun2/analysis2016__SUS-16-028__12p9fb/stopBabies__v8.0.x_v8__20160729/JESdn/";


	if( runlooper || runlostlep || run1lw ) {

		// Data samples
		data->AddFile( dataPath + "data_met_Run2016*_MINIAOD_PromptReco-v2*.root" );
		data->AddFile( dataPath + "data_single_electron_Run2016*_MINIAOD_PromptReco-v2*.root" );
		data->AddFile( dataPath + "data_single_muon_Run2016*_MINIAOD_PromptReco-v2*.root" );

		// Signal sample(s)
		signal->AddFile( sigPath + "Signal_T2tt*.root" );

		// Background samples
		tt2l->AddFile( bkgPath + "ttbar_diLept_madgraph_pythia8_ext1_25ns*.root" );

		tt1l->AddFile( bkgPath + "ttbar_singleLeptFromT_madgraph_pythia8_*.root" );
		tt1l->AddFile( bkgPath + "ttbar_singleLeptFromTbar_madgraph_pythia8_ext1*.root" );

		wjets->AddFile( bkgPath + "WJetsToLNu_HT100To200_madgraph_pythia8_ext1_25ns*.root" );
		wjets->AddFile( bkgPath + "WJetsToLNu_HT200To400_madgraph_pythia8_ext1_25ns*.root" );
		wjets->AddFile( bkgPath + "WJetsToLNu_HT400To600_madgraph_pythia8_25ns.root" ); // extended sample available
		wjets->AddFile( bkgPath + "WJetsToLNu_HT600To800_madgraph_pythia8_25ns.root" ); // extended sample available
		wjets->AddFile( bkgPath + "WJetsToLNu_HT800To1200_madgraph_pythia8_ext1_25ns*.root" );
		wjets->AddFile( bkgPath + "WJetsToLNu_HT1200To2500_madgraph_pythia8_25ns.root" );
		wjets->AddFile( bkgPath + "WJetsToLNu_HT2500ToInf_madgraph_pythia8_25ns.root" ); // extended sample available

		dy->AddFile( bkgPath + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns.root" );
		dy->AddFile( bkgPath + "DYJetsToLL_m50_amcnlo_pythia8_25ns.root" );

		singtop->AddFile( bkgPath + "t_sch_4f_amcnlo_pythia8_25ns.root" );
		// singtop->AddFile( bkgPath + "t_tch_4f_powheg_pythia8_25ns.root" );
		// singtop->AddFile( bkgPath + "tbar_tch_4f_powheg_pythia8_25ns.root" );
		singtop->AddFile( bkgPath + "t_tW_5f_powheg_pythia8_25ns.root" );
		singtop->AddFile( bkgPath + "t_tbarW_5f_powheg_pythia8_25ns.root" );

		rare->AddFile( bkgPath + "TTWJetsToLNu_amcnlo_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "TTWJetsToQQ_amcnlo_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "ttZJets_13TeV_madgraphMLM*.root" );
		// rare->AddFile( bkgPath + "tZq_ll_4f_amcnlo_pythia8_25ns.root" );
		// rare->AddFile( bkgPath + "tZq_nunu_4f_amcnlo_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "WWTo2l2Nu_powheg_25ns.root" );
		rare->AddFile( bkgPath + "WWToLNuQQ_powheg_25ns.root" );
		rare->AddFile( bkgPath + "WZTo3LNu_powheg_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "WZTo2L2Q_amcnlo_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "WZTo1LNu2Q_amcnlo_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "WZTo1L3Nu_amcnlo_pythia8_25ns.root" );
		// rare->AddFile( bkgPath + "ZZTo4L_powheg_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "ZZTo2L2Q_amcnlo_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "ZZTo2L2Nu_powheg_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "ZZTo2Q2Nu_amcnlo_pythia8_25ns.root" );
	}

	if( runjes ) {
		signal_jesup->AddFile( sigPath_jesup + "Signal_T2tt*.root" );
		tt2l_jesup->AddFile( bkgPath_jesup + "ttbar_diLept_madgraph_pythia8_ext1_25ns*.root" );
		tt1l_jesup->AddFile( bkgPath_jesup + "ttbar_singleLeptFromT_madgraph_pythia8_*.root" );
		tt1l_jesup->AddFile( bkgPath_jesup + "ttbar_singleLeptFromTbar_madgraph_pythia8_ext1*.root" );
		wjets_jesup->AddFile( bkgPath_jesup + "WJetsToLNu_HT100To200_madgraph_pythia8_ext1_25ns*.root" );
		wjets_jesup->AddFile( bkgPath_jesup + "WJetsToLNu_HT200To400_madgraph_pythia8_ext1_25ns*.root" );
		wjets_jesup->AddFile( bkgPath_jesup + "WJetsToLNu_HT400To600_madgraph_pythia8_25ns.root" ); // extended sample available
		wjets_jesup->AddFile( bkgPath_jesup + "WJetsToLNu_HT600To800_madgraph_pythia8_25ns.root" ); // extended sample available
		wjets_jesup->AddFile( bkgPath_jesup + "WJetsToLNu_HT800To1200_madgraph_pythia8_ext1_25ns*.root" );
		wjets_jesup->AddFile( bkgPath_jesup + "WJetsToLNu_HT1200To2500_madgraph_pythia8_25ns.root" );
		wjets_jesup->AddFile( bkgPath_jesup + "WJetsToLNu_HT2500ToInf_madgraph_pythia8_25ns.root" ); // extended sample available
		dy_jesup->AddFile( bkgPath_jesup + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns.root" );
		dy_jesup->AddFile( bkgPath_jesup + "DYJetsToLL_m50_amcnlo_pythia8_25ns.root" );
		singtop_jesup->AddFile( bkgPath_jesup + "t_sch_4f_amcnlo_pythia8_25ns.root" );
		// singtop_jesup->AddFile( bkgPath_jesup + "t_tch_4f_powheg_pythia8_25ns.root" );
		// singtop_jesup->AddFile( bkgPath_jesup + "tbar_tch_4f_powheg_pythia8_25ns.root" );
		singtop_jesup->AddFile( bkgPath_jesup + "t_tW_5f_powheg_pythia8_25ns.root" );
		singtop_jesup->AddFile( bkgPath_jesup + "t_tbarW_5f_powheg_pythia8_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "TTWJetsToLNu_amcnlo_pythia8_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "TTWJetsToQQ_amcnlo_pythia8_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "ttZJets_13TeV_madgraphMLM*.root" );
		// rare_jesup->AddFile( bkgPath_jesup + "tZq_ll_4f_amcnlo_pythia8_25ns.root" );
		// rare_jesup->AddFile( bkgPath_jesup + "tZq_nunu_4f_amcnlo_pythia8_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "WWTo2l2Nu_powheg_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "WWToLNuQQ_powheg_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "WZTo3LNu_powheg_pythia8_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "WZTo2L2Q_amcnlo_pythia8_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "WZTo1LNu2Q_amcnlo_pythia8_25ns*.root" );
		rare_jesup->AddFile( bkgPath_jesup + "WZTo1L3Nu_amcnlo_pythia8_25ns.root" );
		// rare_jesup->AddFile( bkgPath_jesup + "ZZTo4L_powheg_pythia8_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "ZZTo2L2Q_amcnlo_pythia8_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "ZZTo2L2Nu_powheg_pythia8_25ns.root" );
		rare_jesup->AddFile( bkgPath_jesup + "ZZTo2Q2Nu_amcnlo_pythia8_25ns.root" );

		signal_jesdn->AddFile( sigPath_jesdn + "Signal_T2tt*.root" );
		tt2l_jesdn->AddFile( bkgPath_jesdn + "ttbar_diLept_madgraph_pythia8_ext1_25ns*.root" );
		tt1l_jesdn->AddFile( bkgPath_jesdn + "ttbar_singleLeptFromT_madgraph_pythia8_*.root" );
		tt1l_jesdn->AddFile( bkgPath_jesdn + "ttbar_singleLeptFromTbar_madgraph_pythia8_ext1*.root" );
		wjets_jesdn->AddFile( bkgPath_jesdn + "WJetsToLNu_HT100To200_madgraph_pythia8_ext1_25ns*.root" );
		wjets_jesdn->AddFile( bkgPath_jesdn + "WJetsToLNu_HT200To400_madgraph_pythia8_ext1_25ns*.root" );
		wjets_jesdn->AddFile( bkgPath_jesdn + "WJetsToLNu_HT400To600_madgraph_pythia8_25ns.root" ); // extended sample available
		wjets_jesdn->AddFile( bkgPath_jesdn + "WJetsToLNu_HT600To800_madgraph_pythia8_25ns.root" ); // extended sample available
		wjets_jesdn->AddFile( bkgPath_jesdn + "WJetsToLNu_HT800To1200_madgraph_pythia8_ext1_25ns*.root" );
		wjets_jesdn->AddFile( bkgPath_jesdn + "WJetsToLNu_HT1200To2500_madgraph_pythia8_25ns.root" );
		wjets_jesdn->AddFile( bkgPath_jesdn + "WJetsToLNu_HT2500ToInf_madgraph_pythia8_25ns.root" ); // extended sample available
		dy_jesdn->AddFile( bkgPath_jesdn + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns.root" );
		dy_jesdn->AddFile( bkgPath_jesdn + "DYJetsToLL_m50_amcnlo_pythia8_25ns.root" );
		singtop_jesdn->AddFile( bkgPath_jesdn + "t_sch_4f_amcnlo_pythia8_25ns.root" );
		// singtop_jesdn->AddFile( bkgPath_jesdn + "t_tch_4f_powheg_pythia8_25ns.root" );
		// singtop_jesdn->AddFile( bkgPath_jesdn + "tbar_tch_4f_powheg_pythia8_25ns.root" );
		singtop_jesdn->AddFile( bkgPath_jesdn + "t_tW_5f_powheg_pythia8_25ns.root" );
		singtop_jesdn->AddFile( bkgPath_jesdn + "t_tbarW_5f_powheg_pythia8_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "TTWJetsToLNu_amcnlo_pythia8_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "TTWJetsToQQ_amcnlo_pythia8_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "ttZJets_13TeV_madgraphMLM*.root" );
		// rare_jesdn->AddFile( bkgPath_jesdn + "tZq_ll_4f_amcnlo_pythia8_25ns.root" );
		// rare_jesdn->AddFile( bkgPath_jesdn + "tZq_nunu_4f_amcnlo_pythia8_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "WWTo2l2Nu_powheg_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "WWToLNuQQ_powheg_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "WZTo3LNu_powheg_pythia8_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "WZTo2L2Q_amcnlo_pythia8_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "WZTo1LNu2Q_amcnlo_pythia8_25ns*.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "WZTo1L3Nu_amcnlo_pythia8_25ns.root" );
		// rare_jesdn->AddFile( bkgPath_jesdn + "ZZTo4L_powheg_pythia8_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "ZZTo2L2Q_amcnlo_pythia8_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "ZZTo2L2Nu_powheg_pythia8_25ns.root" );
		rare_jesdn->AddFile( bkgPath_jesdn + "ZZTo2Q2Nu_amcnlo_pythia8_25ns.root" );
	}

	// Make the 0-bjet control region "analyses" by copying the signal region "analyses" and changing a few properties
	cr0bjets = srAnalysis->Copy( 12.9, "plots0b.root", "systVariations0b.root" );
	cr0b_jesup = sr_jesup->Copy( 12.9, "jes_cr0b.root", "jes_cr0b.root" );
	cr0b_jesdn = sr_jesdn->Copy( 12.9, "jes_cr0b.root", "jes_cr0b.root" );
	cr0bjets->AddSample( data );
	cr0bjets->ResetSystematics();
	cr0bjets->AddSystematics( {&jesup, &jesdn, &lepSFup, &lepSFdn, &btagHFup, &btagHFdn, &btagLFup, &btagLFdn} );
	cr0bjets->AddSystematics( {&qSquaredup, &qSquareddn, &alphaSup, &alphaSdn, &metresup, &metresdn, &contam1lwup, &contam1lwdn } );



	////////////////////////////////////////////////
	// Run ScanChain (the signal region looper)

	if( runlooper ) {

		// Reset output file
		TFile* outfile = new TFile( srAnalysis->GetPlotFileName(), "RECREATE" );
		outfile->Close();
		outfile = new TFile( srAnalysis->GetSystFileName(), "RECREATE" );
		outfile->Close();

		// Run ScanChain on all samples
		for( sample* mySample : srAnalysis->GetAllSamples() ) ScanChain( srAnalysis, mySample );
	}

	////////////////////////////////////////////////
	// Run looperCR2lep (the lost-lepton CR looper)

	if( runlostlep ) {

		// Reset output file
		TFile* outfile = new TFile( crLostLep->GetPlotFileName(), "RECREATE" );
		outfile->Close();
		outfile = new TFile( crLostLep->GetSystFileName(), "RECREATE" );
		outfile->Close();

		// Run lost lepton CR looper on all samples
		for( sample* mySample : crLostLep->GetAllSamples() ) looperCR2lep( crLostLep, mySample );
	}

	////////////////////////////////////////////////
	// Run looperCR0b (the 0-bjet CR looper)

	if( run1lw ) {

		// Reset output file
		TFile* outfile = new TFile( cr0bjets->GetPlotFileName(), "RECREATE" );
		outfile->Close();
		outfile = new TFile( cr0bjets->GetSystFileName(), "RECREATE" );
		outfile->Close();

		// Run 0-bjet CR looper on all samples
		for( sample* mySample : cr0bjets->GetAllSamples() ) looperCR0b( cr0bjets, mySample );
	}

	////////////////////////////////////////////////
	// Run jet energy scale (JES) systematics

	if( runjes ) {

		// Reset output files
		TFile* outfile = new TFile( sr_jesup->GetPlotFileName(), "RECREATE");
		outfile->Close();
		outfile = new TFile( cr2l_jesup->GetPlotFileName(), "RECREATE");
		outfile->Close();
		outfile = new TFile( cr0b_jesup->GetPlotFileName(), "RECREATE");
		outfile->Close();

		// Run ScanChain, looperCR2lep, and loopercr0b on JES up/down babies
		for( sample* mySample : sr_jesup->GetAllSamples()   ) ScanChain(    sr_jesup,   mySample );
		for( sample* mySample : cr2l_jesup->GetAllSamples() ) looperCR2lep( cr2l_jesup, mySample );
		for( sample* mySample : cr0b_jesup->GetAllSamples() ) looperCR0b(   cr0b_jesup, mySample );
		for( sample* mySample : sr_jesdn->GetAllSamples()   ) ScanChain(    sr_jesdn,   mySample );
		for( sample* mySample : cr2l_jesdn->GetAllSamples() ) looperCR2lep( cr2l_jesdn, mySample );
		for( sample* mySample : cr0b_jesdn->GetAllSamples() ) looperCR0b(   cr0b_jesdn, mySample );
	}


	/////////////////////////////////////////////////////////////////////
	// Make all the various outputs - tables, plots, datacards...

	if( runtables )               makeTables( srAnalysis );
	if( runlostlep || runtables ) makeTables( crLostLep );
	if( run1lw || runtables )     makeTables( cr0bjets );
	if( runstacks )               makeStack( srAnalysis );
	if( runestimate ) {
		makeLostLepEstimate( srAnalysis, crLostLep );
		make1lWEstimate( srAnalysis, cr0bjets );
	}
	if( runcards )                makeDataCards( srAnalysis, crLostLep, cr0bjets );

	// Clean up /////////
	delete srAnalysis;
	delete crLostLep;
	delete cr0bjets;
	delete data;
	delete signal;
	delete tt2l;
	delete tt1l;
	delete singtop;
	delete wjets;
	delete dy;
	delete rare;
	delete signal_jesup;
	delete tt2l_jesup;
	delete tt1l_jesup;
	delete singtop_jesup;
	delete wjets_jesup;
	delete dy_jesup;
	delete rare_jesup;	
	delete signal_jesdn;
	delete tt2l_jesdn;
	delete tt1l_jesdn;
	delete singtop_jesdn;
	delete wjets_jesdn;
	delete dy_jesdn;
	delete rare_jesdn;

	return 0;
}
