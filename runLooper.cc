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
extern double dphilmet;
extern double lep1pt;
extern double myMlb;
extern int nTightTags;

bool j1_isBtag;
double j1pt;
double dphilmet;
double lep1pt;
double myMlb;
int nTightTags;


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
	analysis* srAnalysis = new analysis( 29.53, "plots.root", "systVariations.root" );
	analysis* crLostLep  = new analysis( 29.53, "plotsLL.root", "systVariationsLL.root" );
	analysis* cr0bjets   = new analysis( 29.53, "plots0b.root", "systVariations0b.root" );
	analysis* sr_jesup   = new analysis( 29.53, "jes_sr.root", "jes_sr.root" );
	analysis* sr_jesdn   = new analysis( 29.53, "jes_sr.root", "jes_sr.root" );
	analysis* cr2l_jesup = new analysis( 29.53, "jes_cr2l.root", "jes_cr2l.root" );
	analysis* cr2l_jesdn = new analysis( 29.53, "jes_cr2l.root", "jes_cr2l.root" );
	analysis* cr0b_jesup = new analysis( 29.53, "jes_cr0b.root", "jes_cr0b.root" );
	analysis* cr0b_jesdn = new analysis( 29.53, "jes_cr0b.root", "jes_cr0b.root" );


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
	cr0bjets->AddSample(  data );
	srAnalysis->AddSamples( {signal, tt2l, tt1l, singtop, wjets, dy, rare} );
	crLostLep->AddSamples(  {signal, tt2l, tt1l, singtop, wjets, dy, rare} );
	cr0bjets->AddSamples(   {signal, tt2l, tt1l, singtop, wjets, dy, rare} );

	sr_jesup->AddSamples(   {signal_jesup, tt2l_jesup, tt1l_jesup, singtop_jesup, wjets_jesup, dy_jesup, rare_jesup} );
	cr2l_jesup->AddSamples( {signal_jesup, tt2l_jesup, tt1l_jesup, singtop_jesup, wjets_jesup, dy_jesup, rare_jesup} );
	cr0b_jesup->AddSamples( {signal_jesup, tt2l_jesup, tt1l_jesup, singtop_jesup, wjets_jesup, dy_jesup, rare_jesup} );
	sr_jesdn->AddSamples(   {signal_jesdn, tt2l_jesdn, tt1l_jesdn, singtop_jesdn, wjets_jesdn, dy_jesdn, rare_jesdn} );
	cr2l_jesdn->AddSamples( {signal_jesdn, tt2l_jesdn, tt1l_jesdn, singtop_jesdn, wjets_jesdn, dy_jesdn, rare_jesdn} );
	cr0b_jesdn->AddSamples( {signal_jesdn, tt2l_jesdn, tt1l_jesdn, singtop_jesdn, wjets_jesdn, dy_jesdn, rare_jesdn} );

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
	cr0bjets->AddSystematics( {&jesup, &jesdn, &lepSFup, &lepSFdn, &btagHFup, &btagHFdn, &btagLFup, &btagLFdn} );
	cr0bjets->AddSystematics( {&qSquaredup, &qSquareddn, &alphaSup, &alphaSdn, &metresup, &metresdn, &contam1lwup, &contam1lwdn } );

	// A sneaky trick to make JES systematics work with existing code
	systematic jesup_dummy( "JES", systematic::kUp,   (*sfhelp::Unity) );
	systematic jesdn_dummy( "JES", systematic::kDown, (*sfhelp::Unity) );
	sr_jesup->AddSystematics( {&jesup_dummy} );
	sr_jesdn->AddSystematics( {&jesdn_dummy} );
	cr2l_jesup->AddSystematics( {&jesup_dummy} );
	cr2l_jesdn->AddSystematics( {&jesdn_dummy} );
	cr0b_jesup->AddSystematics( {&jesup_dummy} );
	cr0b_jesdn->AddSystematics( {&jesdn_dummy} );


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
	selection<float> MET_250_450( (*tas::pfmet), 250., 450. );
	selection<float> MET_350_550( (*tas::pfmet), 350., 550. );
	selection<float> MET_450_600( (*tas::pfmet), 450., 600. );
	selection<float> MET_600_inf( (*tas::pfmet), 600., 9999999. );

	selection<float> CR2L_MET_250_350( (*tas::pfmet_rl), 250., 350. );
	selection<float> CR2L_MET_350_450( (*tas::pfmet_rl), 350., 450. );
	selection<float> CR2L_MET_450_550( (*tas::pfmet_rl), 450., 550. );
	selection<float> CR2L_MET_550_650( (*tas::pfmet_rl), 550., 650. );
	selection<float> CR2L_MET_450_inf( (*tas::pfmet_rl), 450., 9999999. ); // (MET+lep2) bins for the 2-lep control regions
	selection<float> CR2L_MET_550_inf( (*tas::pfmet_rl), 550., 9999999. );
	selection<float> CR2L_MET_650_inf( (*tas::pfmet_rl), 650., 9999999. );
	selection<float> CR2L_MET_250_450( (*tas::pfmet_rl), 250., 450. );
	selection<float> CR2L_MET_350_550( (*tas::pfmet_rl), 350., 550. );
	selection<float> CR2L_MET_450_600( (*tas::pfmet_rl), 450., 600. );
	selection<float> CR2L_MET_600_inf( (*tas::pfmet_rl), 600., 9999999. );

	selection<int> nJets23(  (*tas::ngoodjets), 2, 3 );
	selection<int> nJetsGe4( (*tas::ngoodjets), 4, 9999999 ); // NJets bins
	selection<int> nJetsGe5( (*tas::ngoodjets), 5, 9999999 );

	selection<int> oneTightB( &nTightTags, 1, 9999999 ); // Various requirements on the number of tight or medium b-tags
	selection<int> noTightBs( &nTightTags, 0 );
	selection<int> noMediumBs( (*tas::ngoodbtags), 0 );

	selection<float> modTopNeg(   (*tas::topnessMod), -9999999., 0. );
	selection<float> modTopLow(   (*tas::topnessMod),  0., 10. );
	selection<float> modTopHigh(  (*tas::topnessMod), 10., 9999999. );   // Modified topness bins
	selection<float> CR2L_modTopNeg(   (*tas::topnessMod_rl), -9999999., 0. );
	selection<float> CR2L_modTopLow(   (*tas::topnessMod_rl),  0., 10. );
	selection<float> CR2L_modTopHigh(  (*tas::topnessMod_rl), 10., 9999999. );

	selection<double> mlbLt175( &myMlb,    0., 175. );
	selection<double> mlbGe175( &myMlb,  175., 9999999. );  // M_lb bins

	selection<double> j1Pt200( &j1pt, 200., 9999999. ); // Special selections for the corridor regions
	selection<bool>   j1NoTag( &j1_isBtag, false );
	selection<double> dPhilepmet(&dphilmet, 0., 1.5 );
	selection<double> lep1ptLt100( &lep1pt, 0., 100. );


	// Create the "sigRegion" objects that will store the definitions of our signal/control regions
	// sigRegion objName(   "label",        "Nice name(s) for plots/tables",             {selections that define the region} )
	sigRegion j23lowmlb250( "j23lowmlb250", "2-3 jets, high tmod, low Mlb, MET250-350",  {&nJets23, &modTopHigh, &mlbLt175, &MET_250_350} ); // A
	sigRegion j23lowmlb350( "j23lowmlb350", "2-3 jets, high tmod, low Mlb, MET350-450",  {&nJets23, &modTopHigh, &mlbLt175, &MET_350_450} );
	sigRegion j23lowmlb450( "j23lowmlb450", "2-3 jets, high tmod, low Mlb, MET450-600",  {&nJets23, &modTopHigh, &mlbLt175, &MET_450_600} );
	sigRegion j23lowmlb600( "j23lowmlb600", "2-3 jets, high tmod, low Mlb, MET600-inf",  {&nJets23, &modTopHigh, &mlbLt175, &MET_600_inf} );
	sigRegion j23himlb250(  "j23himlb250",  "2-3 jets, high tmod, hi Mlb, MET250-450",   {&nJets23, &modTopHigh, &mlbGe175, &MET_250_450, &oneTightB} ); // B
	sigRegion j23himlb450(  "j23himlb450",  "2-3 jets, high tmod, hi Mlb, MET450-600",   {&nJets23, &modTopHigh, &mlbGe175, &MET_450_600, &oneTightB} );
	sigRegion j23himlb600(  "j23himlb600",  "2-3 jets, high tmod, hi Mlb, MET600-inf",   {&nJets23, &modTopHigh, &mlbGe175, &MET_600_inf, &oneTightB} );
	sigRegion j4negtmodlowmlb250( "j4negtmodlowmlb250", "4 jets, neg tmod, low Mlb, MET250-350",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_250_350} ); // C
	sigRegion j4negtmodlowmlb350( "j4negtmodlowmlb350", "4 jets, neg tmod, low Mlb, MET350-450",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_350_450} );
	sigRegion j4negtmodlowmlb450( "j4negtmodlowmlb450", "4 jets, neg tmod, low Mlb, MET450-550",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_450_550} );
	sigRegion j4negtmodlowmlb550( "j4negtmodlowmlb550", "4 jets, neg tmod, low Mlb, MET550-650",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_550_650} );
	sigRegion j4negtmodlowmlb650( "j4negtmodlowmlb650", "4 jets, neg tmod, low Mlb, MET650-inf",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_650_inf} );
	sigRegion j4negtmodhimlb250(  "j4negtmodhimlb250",  "4 jets, neg tmod, hi Mlb, MET250-350",   {&nJetsGe4, &modTopNeg, &mlbGe175, &MET_250_350, &oneTightB} ); // D
	sigRegion j4negtmodhimlb350(  "j4negtmodhimlb350",  "4 jets, neg tmod, hi Mlb, MET350-450",   {&nJetsGe4, &modTopNeg, &mlbGe175, &MET_350_450, &oneTightB} );
	sigRegion j4negtmodhimlb450(  "j4negtmodhimlb450",  "4 jets, neg tmod, hi Mlb, MET450-550",   {&nJetsGe4, &modTopNeg, &mlbGe175, &MET_450_550, &oneTightB} );
	sigRegion j4negtmodhimlb550(  "j4negtmodhimlb550",  "4 jets, neg tmod, hi Mlb, MET550-inf",   {&nJetsGe4, &modTopNeg, &mlbGe175, &MET_550_inf, &oneTightB} );
	sigRegion j4lowtmodlowmlb250( "j4lowtmodlowmlb250", "4 jets, low tmod, low Mlb, MET250-350",  {&nJetsGe4, &modTopLow, &mlbLt175, &MET_250_350} ); // E
	sigRegion j4lowtmodlowmlb350( "j4lowtmodlowmlb350", "4 jets, low tmod, low Mlb, MET350-550",  {&nJetsGe4, &modTopLow, &mlbLt175, &MET_350_550} );
	sigRegion j4lowtmodlowmlb550( "j4lowtmodlowmlb550", "4 jets, low tmod, low Mlb, MET550-inf",  {&nJetsGe4, &modTopLow, &mlbLt175, &MET_550_inf} );
	sigRegion j4lowtmodhimlb250(  "j4lowtmodhimlb250",  "4 jets, low tmod, hi Mlb, MET250-450",   {&nJetsGe4, &modTopLow, &mlbGe175, &MET_250_450, &oneTightB} ); // F
	sigRegion j4lowtmodhimlb450(  "j4lowtmodhimlb450",  "4 jets, low tmod, hi Mlb, MET450-inf",   {&nJetsGe4, &modTopLow, &mlbGe175, &MET_450_inf, &oneTightB} );
	sigRegion j4hitmodlowmlb250( "j4hitmodlowmlb250", "4 jets, high tmod, low Mlb, MET250-350",  {&nJetsGe4, &modTopHigh, &mlbLt175, &MET_250_350} ); // G
	sigRegion j4hitmodlowmlb350( "j4hitmodlowmlb350", "4 jets, high tmod, low Mlb, MET350-450",  {&nJetsGe4, &modTopHigh, &mlbLt175, &MET_350_450} );
	sigRegion j4hitmodlowmlb450( "j4hitmodlowmlb450", "4 jets, high tmod, low Mlb, MET450-600",  {&nJetsGe4, &modTopHigh, &mlbLt175, &MET_450_600} );
	sigRegion j4hitmodlowmlb600( "j4hitmodlowmlb600", "4 jets, high tmod, low Mlb, MET600-inf",  {&nJetsGe4, &modTopHigh, &mlbLt175, &MET_600_inf} );
	sigRegion j4hitmodhimlb250(  "j4hitmodhimlb250",  "4 jets, high tmod, hi Mlb, MET250-450",   {&nJetsGe4, &modTopHigh, &mlbGe175, &MET_250_450, &oneTightB} ); // H
	sigRegion j4hitmodhimlb450(  "j4hitmodhimlb450",  "4 jets, high tmod, hi Mlb, MET450-inf",   {&nJetsGe4, &modTopHigh, &mlbGe175, &MET_450_inf, &oneTightB} );
	sigRegion inclusive(   "inclusive",   "Inclusive" );
	sigRegion corridor250( "corridor250", "Corridor, low MET",  {&nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag} );
	sigRegion corridor350( "corridor350", "Corridor, mid MET",  {&nJetsGe5, &MET_350_450, &j1Pt200, &j1NoTag} );
	sigRegion corridor450( "corridor450", "Corridor, high MET", {&nJetsGe5, &MET_450_inf, &j1Pt200, &j1NoTag} );
	sigRegion corr250combo( "corr250combo", "Corridor 250 combo", {&nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );
	sigRegion corr350combo( "corr350combo", "Corridor 350 combo", {&nJetsGe5, &MET_350_450, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );
	sigRegion corr450combo( "corr450combo", "Corridor 450 combo", {&nJetsGe5, &MET_450_550, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );
	sigRegion corr550combo( "corr550combo", "Corridor 550 combo", {&nJetsGe5, &MET_550_inf, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );
	sigRegion corrAllcombo( "corrAllcombo", "Corridor All combo", {&nJetsGe5,               &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );

	// Equivalents for the diLepton control regions
	sigRegion j23lowmlb250CR2L( "j23lowmlb250CR2L", "CR2L 2-3 jets, high tmod, low Mlb, MET250-350",  {&nJets23, &CR2L_modTopHigh, &mlbLt175, &CR2L_MET_250_350} ); // A
	sigRegion j23lowmlb350CR2L( "j23lowmlb350CR2L", "CR2L 2-3 jets, high tmod, low Mlb, MET350-450",  {&nJets23, &CR2L_modTopHigh, &mlbLt175, &CR2L_MET_350_450} );
	sigRegion j23lowmlb450CR2L( "j23lowmlb450CR2L", "CR2L 2-3 jets, high tmod, low Mlb, MET450-600",  {&nJets23, &CR2L_modTopHigh, &mlbLt175, &CR2L_MET_450_600} );
	sigRegion j23lowmlb600CR2L( "j23lowmlb600CR2L", "CR2L 2-3 jets, high tmod, low Mlb, MET600-inf",  {&nJets23, &CR2L_modTopHigh, &mlbLt175, &CR2L_MET_600_inf} );
	sigRegion j23himlb250CR2L(  "j23himlb250CR2L",  "CR2L 2-3 jets, high tmod, hi Mlb, MET250-450",   {&nJets23, &CR2L_modTopHigh, &mlbGe175, &CR2L_MET_250_450, &oneTightB} ); // B
	sigRegion j23himlb450CR2L(  "j23himlb450CR2L",  "CR2L 2-3 jets, high tmod, hi Mlb, MET450-600",   {&nJets23, &CR2L_modTopHigh, &mlbGe175, &CR2L_MET_450_600, &oneTightB} );
	sigRegion j23himlb600CR2L(  "j23himlb600CR2L",  "CR2L 2-3 jets, high tmod, hi Mlb, MET600-inf",   {&nJets23, &CR2L_modTopHigh, &mlbGe175, &CR2L_MET_600_inf, &oneTightB} );
	sigRegion j4negtmodlowmlb250CR2L( "j4negtmodlowmlb250CR2L", "CR2L 4 jets, neg tmod, low Mlb, MET250-350",  {&nJetsGe4, &CR2L_modTopNeg, &mlbLt175, &CR2L_MET_250_350} ); // C
	sigRegion j4negtmodlowmlb350CR2L( "j4negtmodlowmlb350CR2L", "CR2L 4 jets, neg tmod, low Mlb, MET350-450",  {&nJetsGe4, &CR2L_modTopNeg, &mlbLt175, &CR2L_MET_350_450} );
	sigRegion j4negtmodlowmlb450CR2L( "j4negtmodlowmlb450CR2L", "CR2L 4 jets, neg tmod, low Mlb, MET450-550",  {&nJetsGe4, &CR2L_modTopNeg, &mlbLt175, &CR2L_MET_450_550} );
	sigRegion j4negtmodlowmlb550CR2L( "j4negtmodlowmlb550CR2L", "CR2L 4 jets, neg tmod, low Mlb, MET550-650",  {&nJetsGe4, &CR2L_modTopNeg, &mlbLt175, &CR2L_MET_550_650} );
	sigRegion j4negtmodlowmlb650CR2L( "j4negtmodlowmlb650CR2L", "CR2L 4 jets, neg tmod, low Mlb, MET650-inf",  {&nJetsGe4, &CR2L_modTopNeg, &mlbLt175, &CR2L_MET_650_inf} );
	sigRegion j4negtmodhimlb250CR2L(  "j4negtmodhimlb250CR2L",  "CR2L 4 jets, neg tmod, hi Mlb, MET250-350",   {&nJetsGe4, &CR2L_modTopNeg, &mlbGe175, &CR2L_MET_250_350, &oneTightB} ); // D
	sigRegion j4negtmodhimlb350CR2L(  "j4negtmodhimlb350CR2L",  "CR2L 4 jets, neg tmod, hi Mlb, MET350-450",   {&nJetsGe4, &CR2L_modTopNeg, &mlbGe175, &CR2L_MET_350_450, &oneTightB} );
	sigRegion j4negtmodhimlb450CR2L(  "j4negtmodhimlb450CR2L",  "CR2L 4 jets, neg tmod, hi Mlb, MET450-550",   {&nJetsGe4, &CR2L_modTopNeg, &mlbGe175, &CR2L_MET_450_550, &oneTightB} );
	sigRegion j4negtmodhimlb550CR2L(  "j4negtmodhimlb550CR2L",  "CR2L 4 jets, neg tmod, hi Mlb, MET550-inf",   {&nJetsGe4, &CR2L_modTopNeg, &mlbGe175, &CR2L_MET_550_inf, &oneTightB} );
	sigRegion j4lowtmodlowmlb250CR2L( "j4lowtmodlowmlb250CR2L", "CR2L 4 jets, low tmod, low Mlb, MET250-350",  {&nJetsGe4, &CR2L_modTopLow, &mlbLt175, &CR2L_MET_250_350} ); // E
	sigRegion j4lowtmodlowmlb350CR2L( "j4lowtmodlowmlb350CR2L", "CR2L 4 jets, low tmod, low Mlb, MET350-550",  {&nJetsGe4, &CR2L_modTopLow, &mlbLt175, &CR2L_MET_350_550} );
	sigRegion j4lowtmodlowmlb550CR2L( "j4lowtmodlowmlb550CR2L", "CR2L 4 jets, low tmod, low Mlb, MET550-inf",  {&nJetsGe4, &CR2L_modTopLow, &mlbLt175, &CR2L_MET_550_inf} );
	sigRegion j4lowtmodhimlb250CR2L(  "j4lowtmodhimlb250CR2L",  "CR2L 4 jets, low tmod, hi Mlb, MET250-450",   {&nJetsGe4, &CR2L_modTopLow, &mlbGe175, &CR2L_MET_250_450, &oneTightB} ); // F
	sigRegion j4lowtmodhimlb450CR2L(  "j4lowtmodhimlb450CR2L",  "CR2L 4 jets, low tmod, hi Mlb, MET450-inf",   {&nJetsGe4, &CR2L_modTopLow, &mlbGe175, &CR2L_MET_450_inf, &oneTightB} );
	sigRegion j4hitmodlowmlb250CR2L( "j4hitmodlowmlb250CR2L", "CR2L 4 jets, high tmod, low Mlb, MET250-350",  {&nJetsGe4, &CR2L_modTopHigh, &mlbLt175, &CR2L_MET_250_350} ); // G
	sigRegion j4hitmodlowmlb350CR2L( "j4hitmodlowmlb350CR2L", "CR2L 4 jets, high tmod, low Mlb, MET350-450",  {&nJetsGe4, &CR2L_modTopHigh, &mlbLt175, &CR2L_MET_350_450} );
	sigRegion j4hitmodlowmlb450CR2L( "j4hitmodlowmlb450CR2L", "CR2L 4 jets, high tmod, low Mlb, MET450-600",  {&nJetsGe4, &CR2L_modTopHigh, &mlbLt175, &CR2L_MET_450_600} );
	sigRegion j4hitmodlowmlb600CR2L( "j4hitmodlowmlb600CR2L", "CR2L 4 jets, high tmod, low Mlb, MET600-inf",  {&nJetsGe4, &CR2L_modTopHigh, &mlbLt175, &CR2L_MET_600_inf} );
	sigRegion j4hitmodhimlb250CR2L(  "j4hitmodhimlb250CR2L",  "CR2L 4 jets, high tmod, hi Mlb, MET250-450",   {&nJetsGe4, &CR2L_modTopHigh, &mlbGe175, &CR2L_MET_250_450, &oneTightB} ); // H
	sigRegion j4hitmodhimlb450CR2L(  "j4hitmodhimlb450CR2L",  "CR2L 4 jets, high tmod, hi Mlb, MET450-inf",   {&nJetsGe4, &CR2L_modTopHigh, &mlbGe175, &CR2L_MET_450_inf, &oneTightB} );
	// sigRegion inclusive(   "inclusive",   "Inclusive" );
	sigRegion corridor250CR2L( "corridor250CR2L", "CR2L Corridor, low MET",  {&nJetsGe5, &CR2L_MET_250_350, &j1Pt200, &j1NoTag} );
	sigRegion corridor350CR2L( "corridor350CR2L", "CR2L Corridor, mid MET",  {&nJetsGe5, &CR2L_MET_350_450, &j1Pt200, &j1NoTag} );
	sigRegion corridor450CR2L( "corridor450CR2L", "CR2L Corridor, high MET", {&nJetsGe5, &CR2L_MET_450_inf, &j1Pt200, &j1NoTag} );
	sigRegion corr250comboCR2L( "corr250comboCR2L", "CR2L Corridor 250 combo", {&nJetsGe5, &CR2L_MET_250_350, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );
	sigRegion corr350comboCR2L( "corr350comboCR2L", "CR2L Corridor 350 combo", {&nJetsGe5, &CR2L_MET_350_450, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );
	sigRegion corr450comboCR2L( "corr450comboCR2L", "CR2L Corridor 450 combo", {&nJetsGe5, &CR2L_MET_450_550, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );
	sigRegion corr550comboCR2L( "corr550comboCR2L", "CR2L Corridor 550 combo", {&nJetsGe5, &CR2L_MET_550_inf, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );
	// sigRegion corrAllcomboCR2L( "corrAllcomboCR2L", "CR2L Corridor All combo", {&nJetsGe5,                  &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet} );

	// Equivalents for the 0-btag control regions
	sigRegion j23lowmlb250CR0b( "j23lowmlb250CR0b", "CR0b 2-3 jets, high tmod, low Mlb, MET250-350",  {&nJets23, &modTopHigh, &mlbLt175, &MET_250_350, &noMediumBs} ); // A
	sigRegion j23lowmlb350CR0b( "j23lowmlb350CR0b", "CR0b 2-3 jets, high tmod, low Mlb, MET350-450",  {&nJets23, &modTopHigh, &mlbLt175, &MET_350_450, &noMediumBs} );
	sigRegion j23lowmlb450CR0b( "j23lowmlb450CR0b", "CR0b 2-3 jets, high tmod, low Mlb, MET450-600",  {&nJets23, &modTopHigh, &mlbLt175, &MET_450_600, &noMediumBs} );
	sigRegion j23lowmlb600CR0b( "j23lowmlb600CR0b", "CR0b 2-3 jets, high tmod, low Mlb, MET600-inf",  {&nJets23, &modTopHigh, &mlbLt175, &MET_600_inf, &noMediumBs} );
	sigRegion j23himlb250CR0b(  "j23himlb250CR0b",  "CR0b 2-3 jets, high tmod, hi Mlb, MET250-450",   {&nJets23, &modTopHigh, &mlbGe175, &MET_250_450, &noTightBs} ); // B
	sigRegion j23himlb450CR0b(  "j23himlb450CR0b",  "CR0b 2-3 jets, high tmod, hi Mlb, MET450-600",   {&nJets23, &modTopHigh, &mlbGe175, &MET_450_600, &noTightBs} );
	sigRegion j23himlb600CR0b(  "j23himlb600CR0b",  "CR0b 2-3 jets, high tmod, hi Mlb, MET600-inf",   {&nJets23, &modTopHigh, &mlbGe175, &MET_600_inf, &noTightBs} );
	sigRegion j4negtmodlowmlb250CR0b( "j4negtmodlowmlb250CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET250-350",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_250_350, &noMediumBs} ); // C
	sigRegion j4negtmodlowmlb350CR0b( "j4negtmodlowmlb350CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET350-450",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_350_450, &noMediumBs} );
	sigRegion j4negtmodlowmlb450CR0b( "j4negtmodlowmlb450CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET450-550",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_450_550, &noMediumBs} );
	sigRegion j4negtmodlowmlb550CR0b( "j4negtmodlowmlb550CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET550-650",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_550_650, &noMediumBs} );
	sigRegion j4negtmodlowmlb650CR0b( "j4negtmodlowmlb650CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET650-inf",  {&nJetsGe4, &modTopNeg, &mlbLt175, &MET_650_inf, &noMediumBs} );
	sigRegion j4negtmodhimlb250CR0b(  "j4negtmodhimlb250CR0b",  "CR0b 4 jets, neg tmod, hi Mlb, MET250-350",   {&nJetsGe4, &modTopNeg, &mlbGe175, &MET_250_350, &noTightBs} ); // D
	sigRegion j4negtmodhimlb350CR0b(  "j4negtmodhimlb350CR0b",  "CR0b 4 jets, neg tmod, hi Mlb, MET350-450",   {&nJetsGe4, &modTopNeg, &mlbGe175, &MET_350_450, &noTightBs} );
	sigRegion j4negtmodhimlb450CR0b(  "j4negtmodhimlb450CR0b",  "CR0b 4 jets, neg tmod, hi Mlb, MET450-550",   {&nJetsGe4, &modTopNeg, &mlbGe175, &MET_450_550, &noTightBs} );
	sigRegion j4negtmodhimlb550CR0b(  "j4negtmodhimlb550CR0b",  "CR0b 4 jets, neg tmod, hi Mlb, MET550-inf",   {&nJetsGe4, &modTopNeg, &mlbGe175, &MET_550_inf, &noTightBs} );
	sigRegion j4lowtmodlowmlb250CR0b( "j4lowtmodlowmlb250CR0b", "CR0b 4 jets, low tmod, low Mlb, MET250-350",  {&nJetsGe4, &modTopLow, &mlbLt175, &MET_250_350, &noMediumBs} ); // E
	sigRegion j4lowtmodlowmlb350CR0b( "j4lowtmodlowmlb350CR0b", "CR0b 4 jets, low tmod, low Mlb, MET350-550",  {&nJetsGe4, &modTopLow, &mlbLt175, &MET_350_550, &noMediumBs} );
	sigRegion j4lowtmodlowmlb550CR0b( "j4lowtmodlowmlb550CR0b", "CR0b 4 jets, low tmod, low Mlb, MET550-inf",  {&nJetsGe4, &modTopLow, &mlbLt175, &MET_550_inf, &noMediumBs} );
	sigRegion j4lowtmodhimlb250CR0b(  "j4lowtmodhimlb250CR0b",  "CR0b 4 jets, low tmod, hi Mlb, MET250-450",   {&nJetsGe4, &modTopLow, &mlbGe175, &MET_250_450, &noTightBs} ); // F
	sigRegion j4lowtmodhimlb450CR0b(  "j4lowtmodhimlb450CR0b",  "CR0b 4 jets, low tmod, hi Mlb, MET450-inf",   {&nJetsGe4, &modTopLow, &mlbGe175, &MET_450_inf, &noTightBs} );
	sigRegion j4hitmodlowmlb250CR0b( "j4hitmodlowmlb250CR0b", "CR0b 4 jets, high tmod, low Mlb, MET250-350",  {&nJetsGe4, &modTopHigh, &mlbLt175, &MET_250_350, &noMediumBs} ); // G
	sigRegion j4hitmodlowmlb350CR0b( "j4hitmodlowmlb350CR0b", "CR0b 4 jets, high tmod, low Mlb, MET350-450",  {&nJetsGe4, &modTopHigh, &mlbLt175, &MET_350_450, &noMediumBs} );
	sigRegion j4hitmodlowmlb450CR0b( "j4hitmodlowmlb450CR0b", "CR0b 4 jets, high tmod, low Mlb, MET450-600",  {&nJetsGe4, &modTopHigh, &mlbLt175, &MET_450_600, &noMediumBs} );
	sigRegion j4hitmodlowmlb600CR0b( "j4hitmodlowmlb600CR0b", "CR0b 4 jets, high tmod, low Mlb, MET600-inf",  {&nJetsGe4, &modTopHigh, &mlbLt175, &MET_600_inf, &noMediumBs} );
	sigRegion j4hitmodhimlb250CR0b(  "j4hitmodhimlb250CR0b",  "CR0b 4 jets, high tmod, hi Mlb, MET250-450",   {&nJetsGe4, &modTopHigh, &mlbGe175, &MET_250_450, &noTightBs} ); // H
	sigRegion j4hitmodhimlb450CR0b(  "j4hitmodhimlb450CR0b",  "CR0b 4 jets, high tmod, hi Mlb, MET450-inf",   {&nJetsGe4, &modTopHigh, &mlbGe175, &MET_450_inf, &noTightBs} );
	// sigRegion inclusive(   "inclusive",   "Inclusive" );
	sigRegion corridor250CR0b( "corridor250CR0b", "CR0b Corridor, low MET",  {&nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag, &noMediumBs} );
	sigRegion corridor350CR0b( "corridor350CR0b", "CR0b Corridor, mid MET",  {&nJetsGe5, &MET_350_450, &j1Pt200, &j1NoTag, &noMediumBs} );
	sigRegion corridor450CR0b( "corridor450CR0b", "CR0b Corridor, high MET", {&nJetsGe5, &MET_450_inf, &j1Pt200, &j1NoTag, &noMediumBs} );
	sigRegion corr250comboCR0b( "corr250comboCR0b", "CR0b Corridor 250 combo", {&nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet, &noMediumBs} );
	sigRegion corr350comboCR0b( "corr350comboCR0b", "CR0b Corridor 350 combo", {&nJetsGe5, &MET_350_450, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet, &noMediumBs} );
	sigRegion corr450comboCR0b( "corr450comboCR0b", "CR0b Corridor 450 combo", {&nJetsGe5, &MET_450_550, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet, &noMediumBs} );
	sigRegion corr550comboCR0b( "corr550comboCR0b", "CR0b Corridor 550 combo", {&nJetsGe5, &MET_550_inf, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet, &noMediumBs} );
	sigRegion corrAllcomboCR0b( "corrAllcomboCR0b", "CR0b Corridor All combo", {&nJetsGe5,               &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilepmet, &noMediumBs} );


	// Finally, store all these signal/control regions in our "analysis" objects.
	// Regions grouped together here will also be grouped together in various tables down the road (yields, systematics, etc.)
	srAnalysis->AddSigRegs( {&j23lowmlb250, &j23lowmlb350, &j23lowmlb450, &j23lowmlb600} );
	srAnalysis->AddSigRegs( {&j23himlb250, &j23himlb450, &j23himlb600} );
	srAnalysis->AddSigRegs( {&j4negtmodlowmlb250, &j4negtmodlowmlb350, &j4negtmodlowmlb450, &j4negtmodlowmlb550, &j4negtmodlowmlb650} );
	srAnalysis->AddSigRegs( {&j4negtmodhimlb250, &j4negtmodhimlb350, &j4negtmodhimlb450, &j4negtmodhimlb550} );
	srAnalysis->AddSigRegs( {&j4lowtmodlowmlb250, &j4lowtmodlowmlb350, &j4lowtmodlowmlb550} );
	srAnalysis->AddSigRegs( {&j4lowtmodhimlb250, &j4lowtmodhimlb450} );
	srAnalysis->AddSigRegs( {&j4hitmodlowmlb250, &j4hitmodlowmlb350, &j4hitmodlowmlb450, &j4hitmodlowmlb600} );
	srAnalysis->AddSigRegs( {&j4hitmodhimlb250, &j4hitmodhimlb450} );
	srAnalysis->AddSigRegs( {&inclusive} );
	srAnalysis->AddSigRegs( {&corridor250, &corridor350, &corridor450} );
	srAnalysis->AddSigRegs( {&corr250combo, &corr350combo, &corr450combo, &corr550combo} );
	srAnalysis->AddSigRegs( {&corrAllcombo} );

	crLostLep->AddSigRegs( {&j23lowmlb250CR2L, &j23lowmlb350CR2L, &j23lowmlb450CR2L, &j23lowmlb600CR2L} );
	crLostLep->AddSigRegs( {&j23himlb250CR2L, &j23himlb450CR2L, &j23himlb600CR2L} );
	crLostLep->AddSigRegs( {&j4negtmodlowmlb250CR2L, &j4negtmodlowmlb350CR2L, &j4negtmodlowmlb450CR2L, &j4negtmodlowmlb550CR2L, &j4negtmodlowmlb650CR2L} );
	crLostLep->AddSigRegs( {&j4negtmodhimlb250CR2L, &j4negtmodhimlb350CR2L, &j4negtmodhimlb450CR2L, &j4negtmodhimlb550CR2L} );
	crLostLep->AddSigRegs( {&j4lowtmodlowmlb250CR2L, &j4lowtmodlowmlb350CR2L, &j4lowtmodlowmlb550CR2L} );
	crLostLep->AddSigRegs( {&j4lowtmodhimlb250CR2L, &j4lowtmodhimlb450CR2L} );
	crLostLep->AddSigRegs( {&j4hitmodlowmlb250CR2L, &j4hitmodlowmlb350CR2L, &j4hitmodlowmlb450CR2L, &j4hitmodlowmlb600CR2L} );
	crLostLep->AddSigRegs( {&j4hitmodhimlb250CR2L, &j4hitmodhimlb450CR2L} );
	crLostLep->AddSigRegs( {&inclusive} );
	crLostLep->AddSigRegs( {&corridor250CR2L, &corridor350CR2L, &corridor450CR2L} );
	crLostLep->AddSigRegs( {&corr250comboCR2L, &corr350comboCR2L, &corr450comboCR2L, &corr550comboCR2L} );
	crLostLep->AddSigRegs( {&corrAllcombo} );

	cr0bjets->AddSigRegs( {&j23lowmlb250CR0b, &j23lowmlb350CR0b, &j23lowmlb450CR0b, &j23lowmlb600CR0b} );
	cr0bjets->AddSigRegs( {&j23himlb250CR0b, &j23himlb450CR0b, &j23himlb600CR0b} );
	cr0bjets->AddSigRegs( {&j4negtmodlowmlb250CR0b, &j4negtmodlowmlb350CR0b, &j4negtmodlowmlb450CR0b, &j4negtmodlowmlb550CR0b, &j4negtmodlowmlb650CR0b} );
	cr0bjets->AddSigRegs( {&j4negtmodhimlb250CR0b, &j4negtmodhimlb350CR0b, &j4negtmodhimlb450CR0b, &j4negtmodhimlb550CR0b} );
	cr0bjets->AddSigRegs( {&j4lowtmodlowmlb250CR0b, &j4lowtmodlowmlb350CR0b, &j4lowtmodlowmlb550CR0b} );
	cr0bjets->AddSigRegs( {&j4lowtmodhimlb250CR0b, &j4lowtmodhimlb450CR0b} );
	cr0bjets->AddSigRegs( {&j4hitmodlowmlb250CR0b, &j4hitmodlowmlb350CR0b, &j4hitmodlowmlb450CR0b, &j4hitmodlowmlb600CR0b} );
	cr0bjets->AddSigRegs( {&j4hitmodhimlb250CR0b, &j4hitmodhimlb450CR0b} );
	cr0bjets->AddSigRegs( {&inclusive} );
	cr0bjets->AddSigRegs( {&corridor250CR0b, &corridor350CR0b, &corridor450CR0b} );
	cr0bjets->AddSigRegs( {&corr250comboCR0b, &corr350comboCR0b, &corr450comboCR0b, &corr550comboCR0b} );
	cr0bjets->AddSigRegs( {&corrAllcomboCR0b} );

	// Copy signal regions to JES up and down "analysis" objects
	for( vector<sigRegion*> srList : srAnalysis->GetSigRegions() ) {
		sr_jesup->AddSigRegs( srList );
		sr_jesdn->AddSigRegs( srList );
	}
	for( vector<sigRegion*> srList : crLostLep->GetSigRegions() ) {
		cr2l_jesup->AddSigRegs( srList );
		cr2l_jesdn->AddSigRegs( srList );
	}
	for( vector<sigRegion*> srList : cr0bjets->GetSigRegions() ) {
		cr0b_jesup->AddSigRegs( srList );
		cr0b_jesdn->AddSigRegs( srList );
	}



	//////////////////////////////////////////////////////////////////////////////////////////////
	// For each "sample" object defined earlier, chain up the baby files that make up that sample

	TString sigPath  = "/nfs-7/userdata/stopRun2/analysis2016__SUS-16-028__12p9fb/stopBabies__v8.0.x_v8__20160729/Nominal/";
	TString bkgPath  = "/nfs-7/userdata/isuarez/tupler_babies/merged/Stop_1l/v11/skim/";
	TString dataPath = "/nfs-7/userdata/jgwood/tupler_babies/merged/Stop_1l/v12/skim/";

	TString sigPath_jesup = "/nfs-7/userdata/stopRun2/analysis2016__SUS-16-028__12p9fb/stopBabies__v8.0.x_v8__20160729/JESup/";
	TString bkgPath_jesup = "/nfs-7/userdata/stopRun2/analysis2016__SUS-16-028__12p9fb/stopBabies__v8.0.x_v8__20160729/JESup/";
	TString sigPath_jesdn = "/nfs-7/userdata/stopRun2/analysis2016__SUS-16-028__12p9fb/stopBabies__v8.0.x_v8__20160729/JESdn/";
	TString bkgPath_jesdn = "/nfs-7/userdata/stopRun2/analysis2016__SUS-16-028__12p9fb/stopBabies__v8.0.x_v8__20160729/JESdn/";


	if( runlooper || runlostlep || run1lw ) {

		// Data samples
		data->AddFile( dataPath + "data_met_Run2016*_MINIAOD_PromptReco-v*.root" );
		data->AddFile( dataPath + "data_single_electron_Run2016*_MINIAOD_PromptReco-v*.root" );
		data->AddFile( dataPath + "data_single_muon_Run2016*_MINIAOD_PromptReco-v*.root" );

		// Signal sample(s)
		signal->AddFile( sigPath + "Signal_T2tt*.root" );

		// Background samples
		tt2l->AddFile( bkgPath + "ttbar_diLept_madgraph_pythia8_ext1_25ns*.root" );

		tt1l->AddFile( bkgPath + "ttbar_singleLeptFromT_madgraph_pythia8_ext1*.root" );
		tt1l->AddFile( bkgPath + "ttbar_singleLeptFromTbar_madgraph_pythia8_ext1*.root" );

		wjets->AddFile( bkgPath + "W1JetsToLNu_madgraph_pythia8_25ns.root" );
		wjets->AddFile( bkgPath + "W2JetsToLNu_madgraph_pythia8_25ns.root" ); // Cut out genmet < 200 section from these
		wjets->AddFile( bkgPath + "W3JetsToLNu_madgraph_pythia8_25ns.root" );
		wjets->AddFile( bkgPath + "W4JetsToLNu_madgraph_pythia8_25ns.root" );
		wjets->AddFile( bkgPath + "W1JetsToLNu_nupT200_madgraph_pythia8_25ns.root" );
		wjets->AddFile( bkgPath + "W2JetsToLNu_nupT200_madgraph_pythia8_25ns.root" );
		wjets->AddFile( bkgPath + "W3JetsToLNu_nupT200_madgraph_pythia8_25ns.root" );
		wjets->AddFile( bkgPath + "W4JetsToLNu_nupT200_madgraph_pythia8_25ns.root" );

		dy->AddFile( bkgPath + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns.root" );
		dy->AddFile( bkgPath + "DYJetsToLL_m50_amcnlo_pythia8_25ns.root" );

		singtop->AddFile( bkgPath + "t_sch_4f_amcnlo_pythia8_25ns.root" );
		singtop->AddFile( bkgPath + "t_tch_4f_powheg_pythia8_25ns.root" );
		singtop->AddFile( dataPath + "tbar_tch_4f_powheg_pythia8_25ns*.root" );   // Take this one from John's directory
		singtop->AddFile( dataPath + "t_tW_5f_powheg_pythia8_25ns.root" );       // This one too
		singtop->AddFile( bkgPath + "t_tbarW_5f_powheg_pythia8_25ns.root" );

		rare->AddFile( bkgPath + "ttWJets_13TeV_madgraphMLM.root" );
		rare->AddFile( bkgPath + "ttZJets_13TeV_madgraphMLM*.root" );
		// rare->AddFile( bkgPath + "tZq_ll_4f_amcnlo_pythia8_25ns.root" );
		// rare->AddFile( bkgPath + "tZq_nunu_4f_amcnlo_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "WWTo2l2Nu_powheg_25ns.root" );
		rare->AddFile( bkgPath + "WWToLNuQQ_powheg_25ns.root" );
		rare->AddFile( bkgPath + "WZTo3LNu_powheg_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "WZTo2L2Q_amcnlo_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "WZTo1LNu2Q_amcnlo_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "WZTo1L3Nu_amcnlo_pythia8_25ns.root" );
		rare->AddFile( bkgPath + "ZZTo4L_powheg_pythia8_25ns.root" );
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
