#include "analysis.h"
#include "contextVars.h"
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
extern double lep1pt;
extern double myMlb;

bool j1_isBtag;
double j1pt;
double lep1pt;
double myMlb;

void resetFile( TString filename );


// Help with program options
void printHelp() {
	std::cout << "\nUsage: ./runLooper [arg(s)]\n" << std::endl;

	std::cout << "Takes zero or more arguments (separated by spaces) from the following list:" << std::endl;
	std::cout << "[blank]     equivalent to 'all'" << std::endl;
	std::cout << "help        show this message" << std::endl;
	std::cout << "all         run this entire analysis code" << std::endl;
	std::cout << "scan        run ScanChain only" << std::endl;
	std::cout << "lostlep     run lost lepton looper only" << std::endl;
	std::cout << "1lw         run 1l-from-W background looper only" << std::endl;
	std::cout << "jes         run loopers to calculate JES systematic variations" << std::endl;
	std::cout << "plots       run makeStack only" << std::endl;
	std::cout << "tables      run makeTables only" << std::endl;
	std::cout << "estimate    run background and signal estimates only" << std::endl;
	std::cout << "cards       run makeDataCards only" << std::endl;
	std::cout << "output      run makeStack, makeTables, and makeDataCards only" << std::endl;

	std::cout << "\nUsage example: ./runLooper lostlep tables estimate cards\n" << std::endl;
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

	//                     new analysis( lumi, "histogram storage file", "systematic variation storage file" )
	analysis* srAnalysis = new analysis( 36.46, "plots.root", "systVariations.root" );
	analysis* crLostLep  = new analysis( 36.46, "plotsLL.root", "systVariationsLL.root" );
	analysis* cr0bjets   = new analysis( 36.46, "plots0b.root", "systVariations0b.root" );
	analysis* sr_signal  = new analysis( 36.46, "plotsSig.root", "systVariationsSig.root" );
	analysis* sig_genmet = new analysis( 36.46, "plotsSigGenMet.root", "systVariationsSigGenMet.root" );
	analysis* sr_jesup   = new analysis( 36.46, "jes_sr.root", "jes_sr.root" );
	analysis* sr_jesdn   = new analysis( 36.46, "jes_sr.root", "jes_sr.root" );
	analysis* cr2l_jesup = new analysis( 36.46, "jes_cr2l.root", "jes_cr2l.root" );
	analysis* cr2l_jesdn = new analysis( 36.46, "jes_cr2l.root", "jes_cr2l.root" );
	analysis* cr0b_jesup = new analysis( 36.46, "jes_cr0b.root", "jes_cr0b.root" );
	analysis* cr0b_jesdn = new analysis( 36.46, "jes_cr0b.root", "jes_cr0b.root" );


	//                new sample( "Label",  "Display name(s)", TColor,    sampleType )
	sample* data    = new sample( "data",    "Data",           kBlack,    sample::kData );
	sample* signal  = new sample( "signal",  "T2tt",           kBlue+3,   sample::kSignal );
	sample* sigGenMet = new sample("sigGenMet", "T2tt genmet", kBlue+4,   sample::kSignal );
	sample* tt2l    = new sample( "tt2l", "$t\\bar{t} \\rightarrow 2l$", "t#bar{t} #rightarrow 2l", kCyan-3,   sample::kBackground );
	sample* tt1l    = new sample( "tt1l", "$t\\bar{t} \\rightarrow 1l$", "t#bar{t} #rightarrow 1l", kRed-7,    sample::kBackground );
	sample* singtop = new sample( "singletop", "Single Top",   kGreen-4,  sample::kBackground );
	sample* wjets   = new sample( "wjets",   "W+Jets",         kOrange-2, sample::kBackground );
	sample* dy      = new sample( "dy",      "Drell-Yan",      kRed+2,    sample::kBackground );
	sample* rare    = new sample( "rare",    "Rare",           kMagenta-5,sample::kBackground );

	// srAnalysis->AddSample( data );   // Uncomment this line to unblind
	crLostLep->AddSample( data );
	cr0bjets->AddSample(  data );
	sr_signal->AddSample( signal );
	sig_genmet->AddSample( sigGenMet );
	srAnalysis->AddSamples( {tt2l, tt1l, singtop, wjets, dy, rare} );
	crLostLep->AddSamples(  {signal, tt2l, tt1l, singtop, wjets, dy, rare} );
	cr0bjets->AddSamples(   {signal, tt2l, tt1l, singtop, wjets, dy, rare} );

	sr_jesup->AddSamples(   {signal, tt2l, tt1l, singtop, wjets, dy, rare} );
	cr2l_jesup->AddSamples( {signal, tt2l, tt1l, singtop, wjets, dy, rare} );
	cr0b_jesup->AddSamples( {signal, tt2l, tt1l, singtop, wjets, dy, rare} );
	sr_jesdn->AddSamples(   {signal, tt2l, tt1l, singtop, wjets, dy, rare} );
	cr2l_jesdn->AddSamples( {signal, tt2l, tt1l, singtop, wjets, dy, rare} );
	cr0b_jesdn->AddSamples( {signal, tt2l, tt1l, singtop, wjets, dy, rare} );

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
	systematic metresup(       "METres",   systematic::kUp,    (*sfhelp::MetResUp) );
	systematic metresdn(       "METres",   systematic::kDown,  (*sfhelp::MetResDown) );
	systematic topptup(        "topSysPt", systematic::kUp,    (*sfhelp::TopSystPtUp) );
	systematic topptdn(        "topSysPt", systematic::kDown,  (*sfhelp::TopSystPtDown) );
	systematic contam1lwup(    "contam",   systematic::kUp,    (*sfhelp::Contam1lwUp) );
	systematic contam1lwdn(    "contam",   systematic::kDown,  (*sfhelp::Contam1lwDown) );
	systematic isrnjetsup(     "isrnjets", systematic::kUp,    (*sfhelp::ISRnJetsUp) );
	systematic isrnjetsdn(     "isrnjets", systematic::kDown,  (*sfhelp::ISRnJetsDown) );
	systematic lumi(           "lumi",     systematic::kUp,    (*sfhelp::LumiUp) );
	systematic metavgup(       "METavg",   systematic::kSkipUp,   NULL );
	systematic metavgdn(       "METavg",   systematic::kSkipDown, NULL );

	srAnalysis->AddSystematics( {/*&jesup, &jesdn,*/ &lepSFup, &lepSFdn, /*&btagHFup, &btagHFdn, &btagLFup, &btagLFdn,*/ &qSquaredup, &qSquareddn, &alphaSup, &alphaSdn} );
	srAnalysis->AddSystematics( {&metresup, &metresdn, /*&topptup, &topptdn,*/ &isrnjetsup, &isrnjetsdn} );
	sr_signal->AddSystematics(  {/*&jesup, &jesdn,*/ &lepSFup, &lepSFdn, /*&btagHFup, &btagHFdn, &btagLFup, &btagLFdn,*/ &isrnjetsup, &isrnjetsdn, &lumi, &metavgup, &metavgdn} );
	crLostLep->AddSystematics(  {/*&jesup, &jesdn,*/ &lepSFup, &lepSFdn, /*&btagHFup, &btagHFdn, &btagLFup, &btagLFdn,*/ &qSquaredup, &qSquareddn} );
	crLostLep->AddSystematics(  {&alphaSup, &alphaSdn, &eff2lup, &eff2ldn, &metresup, &metresdn, /*&topptup, &topptdn*/ &isrnjetsup, &isrnjetsdn } );
	cr0bjets->AddSystematics( {/*&jesup, &jesdn,*/ &lepSFup, &lepSFdn, /*&btagHFup, &btagHFdn, &btagLFup, &btagLFdn*/ &qSquaredup, &qSquareddn} );
	cr0bjets->AddSystematics( {&alphaSup, &alphaSdn, &metresup, &metresdn, &contam1lwup, &contam1lwdn, &isrnjetsup, &isrnjetsdn } );

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
	selection<float> MET_250_350( (*context::Met), 250., 350. );
	selection<float> MET_350_450( (*context::Met), 350., 450. );
	selection<float> MET_450_550( (*context::Met), 450., 550. );
	selection<float> MET_550_650( (*context::Met), 550., 650. );
	selection<float> MET_450_inf( (*context::Met), 450., 9999999. ); // MET bins for the signal regions
	selection<float> MET_550_inf( (*context::Met), 550., 9999999. );
	selection<float> MET_650_inf( (*context::Met), 650., 9999999. );
	selection<float> MET_250_450( (*context::Met), 250., 450. );
	selection<float> MET_350_550( (*context::Met), 350., 550. );
	selection<float> MET_450_600( (*context::Met), 450., 600. );
	selection<float> MET_600_inf( (*context::Met), 600., 9999999. );

	selection<int> nJets23(  (*context::ngoodjets), 2, 3 );
	selection<int> nJetsGe4( (*context::ngoodjets), 4, 9999999 ); // NJets bins
	selection<int> nJetsGe5( (*context::ngoodjets), 5, 9999999 );

	selection<int> oneTightB( (*context::ntightbtags), 1, 9999999 ); // Various requirements on the number of tight or medium b-tags
	selection<int> noTightBs( (*context::ntightbtags), 0 );
	selection<int> noMediumBs( (*context::ngoodbtags), 0 );

	selection<float> modTopNeg(   (*context::TopnessMod), -9999999., 0. );
	selection<float> modTopLow(   (*context::TopnessMod),  0., 10. );
	selection<float> modTopHigh(  (*context::TopnessMod), 10., 9999999. );   // Modified topness bins

	selection<double> mlbLt175( &myMlb,    0., 175. );
	selection<double> mlbGe175( &myMlb,  175., 9999999. );  // M_lb bins

	selection<float> dPhilMetLt15( (*context::lep1_dphiMET), 0., 1.5 );   // Special selections for the corridor regions
	selection<float> dPhilMetLt20( (*context::lep1_dphiMET), 0., 2.0 );
	selection<float> minDPhi050( (*context::Mindphi_met_j1_j2), 0.50, 9999999. );
	selection<float> minDPhi080( (*context::Mindphi_met_j1_j2), 0.80, 9999999. );
	selection<double> lep1ptLt100( &lep1pt, 0., 100. );
	selection<double> lep1ptLt150( &lep1pt, 0., 150. );
	selection<double> j1Pt200( &j1pt, 200., 9999999. );
	selection<bool> j1NoTag( &j1_isBtag, false );



	// Create the "sigRegion" objects that will store the definitions of our signal regions
	// sigRegion objName(   "label",        "Nice name(s) for plots/tables",             {selections that define the region} )
	sigRegion j23lowmlb250( "j23lowmlb250", "2-3 jets, high tmod, low Mlb, MET250-350",  {&minDPhi080, &nJets23, &modTopHigh, &mlbLt175, &MET_250_350} ); // A
	sigRegion j23lowmlb350( "j23lowmlb350", "2-3 jets, high tmod, low Mlb, MET350-450",  {&minDPhi080, &nJets23, &modTopHigh, &mlbLt175, &MET_350_450} );
	sigRegion j23lowmlb450( "j23lowmlb450", "2-3 jets, high tmod, low Mlb, MET450-600",  {&minDPhi080, &nJets23, &modTopHigh, &mlbLt175, &MET_450_600} );
	sigRegion j23lowmlb600( "j23lowmlb600", "2-3 jets, high tmod, low Mlb, MET600-inf",  {&minDPhi080, &nJets23, &modTopHigh, &mlbLt175, &MET_600_inf} );
	sigRegion j23himlb250(  "j23himlb250",  "2-3 jets, high tmod, hi Mlb, MET250-450",   {&minDPhi080, &nJets23, &modTopHigh, &mlbGe175, &MET_250_450, &oneTightB} ); // B
	sigRegion j23himlb450(  "j23himlb450",  "2-3 jets, high tmod, hi Mlb, MET450-600",   {&minDPhi080, &nJets23, &modTopHigh, &mlbGe175, &MET_450_600, &oneTightB} );
	sigRegion j23himlb600(  "j23himlb600",  "2-3 jets, high tmod, hi Mlb, MET600-inf",   {&minDPhi080, &nJets23, &modTopHigh, &mlbGe175, &MET_600_inf, &oneTightB} );
	sigRegion j4negtmodlowmlb250( "j4negtmodlowmlb250", "4 jets, neg tmod, low Mlb, MET250-350",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_250_350} ); // C
	sigRegion j4negtmodlowmlb350( "j4negtmodlowmlb350", "4 jets, neg tmod, low Mlb, MET350-450",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_350_450} );
	sigRegion j4negtmodlowmlb450( "j4negtmodlowmlb450", "4 jets, neg tmod, low Mlb, MET450-550",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_450_550} );
	sigRegion j4negtmodlowmlb550( "j4negtmodlowmlb550", "4 jets, neg tmod, low Mlb, MET550-650",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_550_650} );
	sigRegion j4negtmodlowmlb650( "j4negtmodlowmlb650", "4 jets, neg tmod, low Mlb, MET650-inf",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_650_inf} );
	sigRegion j4negtmodhimlb250(  "j4negtmodhimlb250",  "4 jets, neg tmod, hi Mlb, MET250-350",   {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbGe175, &MET_250_350, &oneTightB} ); // D
	sigRegion j4negtmodhimlb350(  "j4negtmodhimlb350",  "4 jets, neg tmod, hi Mlb, MET350-450",   {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbGe175, &MET_350_450, &oneTightB} );
	sigRegion j4negtmodhimlb450(  "j4negtmodhimlb450",  "4 jets, neg tmod, hi Mlb, MET450-550",   {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbGe175, &MET_450_550, &oneTightB} );
	sigRegion j4negtmodhimlb550(  "j4negtmodhimlb550",  "4 jets, neg tmod, hi Mlb, MET550-inf",   {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbGe175, &MET_550_inf, &oneTightB} );
	sigRegion j4lowtmodlowmlb250( "j4lowtmodlowmlb250", "4 jets, low tmod, low Mlb, MET250-350",  {&minDPhi080, &nJetsGe4, &modTopLow, &mlbLt175, &MET_250_350} ); // E
	sigRegion j4lowtmodlowmlb350( "j4lowtmodlowmlb350", "4 jets, low tmod, low Mlb, MET350-550",  {&minDPhi080, &nJetsGe4, &modTopLow, &mlbLt175, &MET_350_550} );
	sigRegion j4lowtmodlowmlb550( "j4lowtmodlowmlb550", "4 jets, low tmod, low Mlb, MET550-inf",  {&minDPhi080, &nJetsGe4, &modTopLow, &mlbLt175, &MET_550_inf} );
	sigRegion j4lowtmodhimlb250(  "j4lowtmodhimlb250",  "4 jets, low tmod, hi Mlb, MET250-450",   {&minDPhi080, &nJetsGe4, &modTopLow, &mlbGe175, &MET_250_450, &oneTightB} ); // F
	sigRegion j4lowtmodhimlb450(  "j4lowtmodhimlb450",  "4 jets, low tmod, hi Mlb, MET450-inf",   {&minDPhi080, &nJetsGe4, &modTopLow, &mlbGe175, &MET_450_inf, &oneTightB} );
	sigRegion j4hitmodlowmlb250( "j4hitmodlowmlb250", "4 jets, high tmod, low Mlb, MET250-350",  {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbLt175, &MET_250_350} ); // G
	sigRegion j4hitmodlowmlb350( "j4hitmodlowmlb350", "4 jets, high tmod, low Mlb, MET350-450",  {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbLt175, &MET_350_450} );
	sigRegion j4hitmodlowmlb450( "j4hitmodlowmlb450", "4 jets, high tmod, low Mlb, MET450-600",  {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbLt175, &MET_450_600} );
	sigRegion j4hitmodlowmlb600( "j4hitmodlowmlb600", "4 jets, high tmod, low Mlb, MET600-inf",  {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbLt175, &MET_600_inf} );
	sigRegion j4hitmodhimlb250(  "j4hitmodhimlb250",  "4 jets, high tmod, hi Mlb, MET250-450",   {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbGe175, &MET_250_450, &oneTightB} ); // H
	sigRegion j4hitmodhimlb450(  "j4hitmodhimlb450",  "4 jets, high tmod, hi Mlb, MET450-inf",   {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbGe175, &MET_450_inf, &oneTightB} );
	sigRegion inclusive(   "inclusive",   "Inclusive" );
	sigRegion corr250ichep( "corr250ichep", "Corridor ICHEP 250",  {&minDPhi080, &nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag} );
	sigRegion corr350ichep( "corr350ichep", "Corridor ICHEP 350",  {&minDPhi080, &nJetsGe5, &MET_350_450, &j1Pt200, &j1NoTag} );
	sigRegion corr450ichep( "corr450ichep", "Corridor ICHEP 450", {&minDPhi080, &nJetsGe5, &MET_450_inf, &j1Pt200, &j1NoTag} );
	sigRegion corr250combo( "corr250combo", "Corridor 250 combo", {&minDPhi080, &nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15} );
	sigRegion corr350combo( "corr350combo", "Corridor 350 combo", {&minDPhi080, &nJetsGe5, &MET_350_450, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15} );
	sigRegion corr450combo( "corr450combo", "Corridor 450 combo", {&minDPhi080, &nJetsGe5, &MET_450_550, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15} );
	sigRegion corr550combo( "corr550combo", "Corridor 550 combo", {&minDPhi080, &nJetsGe5, &MET_550_inf, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15} );
	sigRegion corrAllcombo( "corrAllcombo", "Corridor All combo", {&minDPhi080, &nJetsGe5,               &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15} );
	sigRegion corr250new( "corr250new", "Corridor, MET 250-350", {&minDPhi050, &nJetsGe5, &MET_250_350, &lep1ptLt150, &dPhilMetLt20, &j1NoTag} );
	sigRegion corr350new( "corr350new", "Corridor, MET 350-450", {&minDPhi050, &nJetsGe5, &MET_350_450, &lep1ptLt150, &dPhilMetLt20, &j1NoTag} );
	sigRegion corr450new( "corr450new", "Corridor, MET 450-550", {&minDPhi050, &nJetsGe5, &MET_450_550, &lep1ptLt150, &dPhilMetLt20, &j1NoTag} );
	sigRegion corr550new( "corr550new", "Corridor, MET 550-inf", {&minDPhi050, &nJetsGe5, &MET_550_inf, &lep1ptLt150, &dPhilMetLt20, &j1NoTag} );

	// Store these signal regions in our "analysis" object.
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
	srAnalysis->AddSigRegs( {&corr250ichep, &corr350ichep, &corr450ichep} );
	srAnalysis->AddSigRegs( {&corr250combo, &corr350combo, &corr450combo, &corr550combo} );
	srAnalysis->AddSigRegs( {&corrAllcombo} );
	srAnalysis->AddSigRegs( {&corr250new, &corr350new, &corr450new, &corr550new} );

	// Copy all sigRegions from the SR analysis to the lost lepton and signal analyses.
	// Cuts that depend on e.g. MET_rl will be handled automatically by the "contextVars" class
	for( std::vector<sigRegion*> regList : srAnalysis->GetSigRegions() ) {
		sr_signal->AddSigRegs(  regList );
		sig_genmet->AddSigRegs( regList );
		crLostLep->AddSigRegs(  regList );
	}



	// Also create sigRegion objects for the 0-btag control regions
	sigRegion j23lowmlb250CR0b( "j23lowmlb250CR0b", "CR0b 2-3 jets, high tmod, low Mlb, MET250-350",  {&minDPhi080, &nJets23, &modTopHigh, &mlbLt175, &MET_250_350, &noMediumBs} ); // A
	sigRegion j23lowmlb350CR0b( "j23lowmlb350CR0b", "CR0b 2-3 jets, high tmod, low Mlb, MET350-450",  {&minDPhi080, &nJets23, &modTopHigh, &mlbLt175, &MET_350_450, &noMediumBs} );
	sigRegion j23lowmlb450CR0b( "j23lowmlb450CR0b", "CR0b 2-3 jets, high tmod, low Mlb, MET450-600",  {&minDPhi080, &nJets23, &modTopHigh, &mlbLt175, &MET_450_600, &noMediumBs} );
	sigRegion j23lowmlb600CR0b( "j23lowmlb600CR0b", "CR0b 2-3 jets, high tmod, low Mlb, MET600-inf",  {&minDPhi080, &nJets23, &modTopHigh, &mlbLt175, &MET_600_inf, &noMediumBs} );
	sigRegion j23himlb250CR0b(  "j23himlb250CR0b",  "CR0b 2-3 jets, high tmod, hi Mlb, MET250-450",   {&minDPhi080, &nJets23, &modTopHigh, &mlbGe175, &MET_250_450, &noTightBs} ); // B
	sigRegion j23himlb450CR0b(  "j23himlb450CR0b",  "CR0b 2-3 jets, high tmod, hi Mlb, MET450-600",   {&minDPhi080, &nJets23, &modTopHigh, &mlbGe175, &MET_450_600, &noTightBs} );
	sigRegion j23himlb600CR0b(  "j23himlb600CR0b",  "CR0b 2-3 jets, high tmod, hi Mlb, MET600-inf",   {&minDPhi080, &nJets23, &modTopHigh, &mlbGe175, &MET_600_inf, &noTightBs} );
	sigRegion j4negtmodlowmlb250CR0b( "j4negtmodlowmlb250CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET250-350",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_250_350, &noMediumBs} ); // C
	sigRegion j4negtmodlowmlb350CR0b( "j4negtmodlowmlb350CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET350-450",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_350_450, &noMediumBs} );
	sigRegion j4negtmodlowmlb450CR0b( "j4negtmodlowmlb450CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET450-550",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_450_550, &noMediumBs} );
	sigRegion j4negtmodlowmlb550CR0b( "j4negtmodlowmlb550CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET550-650",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_550_650, &noMediumBs} );
	sigRegion j4negtmodlowmlb650CR0b( "j4negtmodlowmlb650CR0b", "CR0b 4 jets, neg tmod, low Mlb, MET650-inf",  {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbLt175, &MET_650_inf, &noMediumBs} );
	sigRegion j4negtmodhimlb250CR0b(  "j4negtmodhimlb250CR0b",  "CR0b 4 jets, neg tmod, hi Mlb, MET250-350",   {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbGe175, &MET_250_350, &noTightBs} ); // D
	sigRegion j4negtmodhimlb350CR0b(  "j4negtmodhimlb350CR0b",  "CR0b 4 jets, neg tmod, hi Mlb, MET350-450",   {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbGe175, &MET_350_450, &noTightBs} );
	sigRegion j4negtmodhimlb450CR0b(  "j4negtmodhimlb450CR0b",  "CR0b 4 jets, neg tmod, hi Mlb, MET450-550",   {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbGe175, &MET_450_550, &noTightBs} );
	sigRegion j4negtmodhimlb550CR0b(  "j4negtmodhimlb550CR0b",  "CR0b 4 jets, neg tmod, hi Mlb, MET550-inf",   {&minDPhi080, &nJetsGe4, &modTopNeg, &mlbGe175, &MET_550_inf, &noTightBs} );
	sigRegion j4lowtmodlowmlb250CR0b( "j4lowtmodlowmlb250CR0b", "CR0b 4 jets, low tmod, low Mlb, MET250-350",  {&minDPhi080, &nJetsGe4, &modTopLow, &mlbLt175, &MET_250_350, &noMediumBs} ); // E
	sigRegion j4lowtmodlowmlb350CR0b( "j4lowtmodlowmlb350CR0b", "CR0b 4 jets, low tmod, low Mlb, MET350-550",  {&minDPhi080, &nJetsGe4, &modTopLow, &mlbLt175, &MET_350_550, &noMediumBs} );
	sigRegion j4lowtmodlowmlb550CR0b( "j4lowtmodlowmlb550CR0b", "CR0b 4 jets, low tmod, low Mlb, MET550-inf",  {&minDPhi080, &nJetsGe4, &modTopLow, &mlbLt175, &MET_550_inf, &noMediumBs} );
	sigRegion j4lowtmodhimlb250CR0b(  "j4lowtmodhimlb250CR0b",  "CR0b 4 jets, low tmod, hi Mlb, MET250-450",   {&minDPhi080, &nJetsGe4, &modTopLow, &mlbGe175, &MET_250_450, &noTightBs} ); // F
	sigRegion j4lowtmodhimlb450CR0b(  "j4lowtmodhimlb450CR0b",  "CR0b 4 jets, low tmod, hi Mlb, MET450-inf",   {&minDPhi080, &nJetsGe4, &modTopLow, &mlbGe175, &MET_450_inf, &noTightBs} );
	sigRegion j4hitmodlowmlb250CR0b( "j4hitmodlowmlb250CR0b", "CR0b 4 jets, high tmod, low Mlb, MET250-350",  {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbLt175, &MET_250_350, &noMediumBs} ); // G
	sigRegion j4hitmodlowmlb350CR0b( "j4hitmodlowmlb350CR0b", "CR0b 4 jets, high tmod, low Mlb, MET350-450",  {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbLt175, &MET_350_450, &noMediumBs} );
	sigRegion j4hitmodlowmlb450CR0b( "j4hitmodlowmlb450CR0b", "CR0b 4 jets, high tmod, low Mlb, MET450-600",  {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbLt175, &MET_450_600, &noMediumBs} );
	sigRegion j4hitmodlowmlb600CR0b( "j4hitmodlowmlb600CR0b", "CR0b 4 jets, high tmod, low Mlb, MET600-inf",  {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbLt175, &MET_600_inf, &noMediumBs} );
	sigRegion j4hitmodhimlb250CR0b(  "j4hitmodhimlb250CR0b",  "CR0b 4 jets, high tmod, hi Mlb, MET250-450",   {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbGe175, &MET_250_450, &noTightBs} ); // H
	sigRegion j4hitmodhimlb450CR0b(  "j4hitmodhimlb450CR0b",  "CR0b 4 jets, high tmod, hi Mlb, MET450-inf",   {&minDPhi080, &nJetsGe4, &modTopHigh, &mlbGe175, &MET_450_inf, &noTightBs} );
	// sigRegion inclusive(   "inclusive",   "Inclusive" );
	sigRegion corr250ichepCR0b( "corr250ichepCR0b", "CR0b Corridor ICHEP 250",  {&minDPhi080, &nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag, &noMediumBs} );
	sigRegion corr350ichepCR0b( "corr350ichepCR0b", "CR0b Corridor ICHEP 350",  {&minDPhi080, &nJetsGe5, &MET_350_450, &j1Pt200, &j1NoTag, &noMediumBs} );
	sigRegion corr450ichepCR0b( "corr450ichepCR0b", "CR0b Corridor ICHEP 450", {&minDPhi080, &nJetsGe5, &MET_450_inf, &j1Pt200, &j1NoTag, &noMediumBs} );
	sigRegion corr250comboCR0b( "corr250comboCR0b", "CR0b Corridor 250 combo", {&minDPhi080, &nJetsGe5, &MET_250_350, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15, &noMediumBs} );
	sigRegion corr350comboCR0b( "corr350comboCR0b", "CR0b Corridor 350 combo", {&minDPhi080, &nJetsGe5, &MET_350_450, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15, &noMediumBs} );
	sigRegion corr450comboCR0b( "corr450comboCR0b", "CR0b Corridor 450 combo", {&minDPhi080, &nJetsGe5, &MET_450_550, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15, &noMediumBs} );
	sigRegion corr550comboCR0b( "corr550comboCR0b", "CR0b Corridor 550 combo", {&minDPhi080, &nJetsGe5, &MET_550_inf, &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15, &noMediumBs} );
	sigRegion corrAllcomboCR0b( "corrAllcomboCR0b", "CR0b Corridor All combo", {&minDPhi080, &nJetsGe5,               &j1Pt200, &j1NoTag, &lep1ptLt100, &dPhilMetLt15, &noMediumBs} );
	sigRegion corr250newCR0b( "corr250newCR0b", "CR0b Corridor, MET 250-350", {&minDPhi050, &nJetsGe5, &MET_250_350, &lep1ptLt150, &dPhilMetLt20, &noMediumBs, &j1NoTag} );
	sigRegion corr350newCR0b( "corr350newCR0b", "CR0b Corridor, MET 350-450", {&minDPhi050, &nJetsGe5, &MET_350_450, &lep1ptLt150, &dPhilMetLt20, &noMediumBs, &j1NoTag} );
	sigRegion corr450newCR0b( "corr450newCR0b", "CR0b Corridor, MET 450-550", {&minDPhi050, &nJetsGe5, &MET_450_550, &lep1ptLt150, &dPhilMetLt20, &noMediumBs, &j1NoTag} );
	sigRegion corr550newCR0b( "corr550newCR0b", "CR0b Corridor, MET 550-inf", {&minDPhi050, &nJetsGe5, &MET_550_inf, &lep1ptLt150, &dPhilMetLt20, &noMediumBs, &j1NoTag} );

	// Store these sigRegions in the CR0b analysis object
	cr0bjets->AddSigRegs( {&j23lowmlb250CR0b, &j23lowmlb350CR0b, &j23lowmlb450CR0b, &j23lowmlb600CR0b} );
	cr0bjets->AddSigRegs( {&j23himlb250CR0b, &j23himlb450CR0b, &j23himlb600CR0b} );
	cr0bjets->AddSigRegs( {&j4negtmodlowmlb250CR0b, &j4negtmodlowmlb350CR0b, &j4negtmodlowmlb450CR0b, &j4negtmodlowmlb550CR0b, &j4negtmodlowmlb650CR0b} );
	cr0bjets->AddSigRegs( {&j4negtmodhimlb250CR0b, &j4negtmodhimlb350CR0b, &j4negtmodhimlb450CR0b, &j4negtmodhimlb550CR0b} );
	cr0bjets->AddSigRegs( {&j4lowtmodlowmlb250CR0b, &j4lowtmodlowmlb350CR0b, &j4lowtmodlowmlb550CR0b} );
	cr0bjets->AddSigRegs( {&j4lowtmodhimlb250CR0b, &j4lowtmodhimlb450CR0b} );
	cr0bjets->AddSigRegs( {&j4hitmodlowmlb250CR0b, &j4hitmodlowmlb350CR0b, &j4hitmodlowmlb450CR0b, &j4hitmodlowmlb600CR0b} );
	cr0bjets->AddSigRegs( {&j4hitmodhimlb250CR0b, &j4hitmodhimlb450CR0b} );
	cr0bjets->AddSigRegs( {&inclusive} );
	cr0bjets->AddSigRegs( {&corr250ichepCR0b, &corr350ichepCR0b, &corr450ichepCR0b} );
	cr0bjets->AddSigRegs( {&corr250comboCR0b, &corr350comboCR0b, &corr450comboCR0b, &corr550comboCR0b} );
	cr0bjets->AddSigRegs( {&corrAllcomboCR0b} );
	cr0bjets->AddSigRegs( {&corr250newCR0b, &corr350newCR0b, &corr450newCR0b, &corr550newCR0b} );


	// Copy signal/control regions to JES up and down "analysis" objects
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

	TString sigPath  = "/nfs-7/userdata/isuarez/tupler_babies/merged/Stop_1l/v17/skim/";
	TString bkgPath  = "/nfs-7/userdata/isuarez/tupler_babies/merged/Stop_1l/v18/skim/";
	TString dataPath = "/nfs-7/userdata/isuarez/tupler_babies/merged/Stop_1l/v17/skim/";


	if( runlooper || runlostlep || run1lw || runjes ) {

		// Data samples
		data->AddFile( dataPath + "data_met_Run2016*_MINIAOD_*.root" );
		data->AddFile( dataPath + "data_single_electron_Run2016*_MINIAOD_*.root" );
		data->AddFile( dataPath + "data_single_muon_Run2016*_MINIAOD_*.root" );

		// Signal sample(s)
		signal->AddFile( sigPath + "Signal_T2tt*.root" );
		sigGenMet->AddFile( "/nfs-7/userdata/haweber/tupler_babies/merged/Stop_1l/v16_tempSignal_gen_oldJEC/output/Signal_T2tt*.root" );

		// Background samples
		tt2l->AddFile( bkgPath + "ttbar_diLept_madgraph_pythia8_ext1_25ns*.root" );

		tt1l->AddFile( bkgPath + "ttbar_singleLeptFromT_madgraph_pythia8_ext1*.root" );
		tt1l->AddFile( bkgPath + "ttbar_singleLeptFromTbar_madgraph_pythia8_ext1*.root" );

		wjets->AddFile( bkgPath + "W1JetsToLNu_madgraph_pythia8_25ns*.root" );
		wjets->AddFile( bkgPath + "W2JetsToLNu_madgraph_pythia8_25ns*.root" ); // Cut out genmet < 200 section from these
		wjets->AddFile( bkgPath + "W3JetsToLNu_madgraph_pythia8_25ns*.root" );
		wjets->AddFile( bkgPath + "W4JetsToLNu_madgraph_pythia8_25ns*.root" );
		wjets->AddFile( bkgPath + "W1JetsToLNu_nupT200_madgraph_pythia8_25ns*.root" );
		wjets->AddFile( bkgPath + "W2JetsToLNu_nupT200_madgraph_pythia8_25ns*.root" );
		wjets->AddFile( bkgPath + "W3JetsToLNu_nupT200_madgraph_pythia8_25ns*.root" );
		wjets->AddFile( bkgPath + "W4JetsToLNu_nupT200_madgraph_pythia8_25ns*.root" );

		dy->AddFile( bkgPath + "DYJetsToLL_m10To50_amcnlo_pythia8_25ns*.root" );
		dy->AddFile( bkgPath + "DYJetsToLL_m50_amcnlo_pythia8_25ns*.root" );

		singtop->AddFile( bkgPath + "t_sch_4f_amcnlo_pythia8_25ns*.root" );
		singtop->AddFile( bkgPath + "t_tch_4f_powheg_pythia8_25ns*.root" );
		singtop->AddFile( bkgPath + "tbar_tch_4f_powheg_pythia8_25ns*.root" );
		singtop->AddFile( bkgPath + "t_tW_5f_powheg_pythia8_25ns*.root" );
		singtop->AddFile( bkgPath + "t_tbarW_5f_powheg_pythia8_25ns*.root" );

		rare->AddFile( bkgPath + "ttWJets_13TeV_madgraphMLM*.root" );
		rare->AddFile( bkgPath + "ttZJets_13TeV_madgraphMLM*.root" );
		// rare->AddFile( bkgPath + "tZq_ll_4f_amcnlo_pythia8_25ns*.root" );
		// rare->AddFile( bkgPath + "tZq_nunu_4f_amcnlo_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "WWTo2l2Nu_powheg_25ns*.root" );
		rare->AddFile( bkgPath + "WWToLNuQQ_powheg_25ns*.root" );
		rare->AddFile( bkgPath + "WZTo3LNu_powheg_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "WZTo2L2Q_amcnlo_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "WZTo1LNu2Q_amcnlo_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "WZTo1L3Nu_amcnlo_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "ZZTo4L_powheg_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "ZZTo2L2Q_amcnlo_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "ZZTo2L2Nu_powheg_pythia8_25ns*.root" );
		rare->AddFile( bkgPath + "ZZTo2Q2Nu_amcnlo_pythia8_25ns*.root" );
	}



	////////////////////////////////////////////////
	// Run ScanChain (the signal region looper)

	if( runlooper ) {

		myContext.SetUseRl( false );
		myContext.SetJesDir( contextVars::kNominal );

		// Run ScanChain on signal sample...
		resetFile( sr_signal->GetPlotFileName() );
		resetFile( sr_signal->GetSystFileName() );
		ScanChain( sr_signal, signal );

		// ...on genmet signal sample...
		resetFile( sig_genmet->GetPlotFileName() );
		resetFile( sig_genmet->GetSystFileName() );
		ScanChain( sig_genmet, sigGenMet );

		// ...and on all other samples
		resetFile( srAnalysis->GetPlotFileName() );
		resetFile( srAnalysis->GetSystFileName() );
		for( sample* mySample : srAnalysis->GetAllSamples() ) ScanChain( srAnalysis, mySample );
	}

	////////////////////////////////////////////////
	// Run looperCR2lep (the lost-lepton CR looper)

	if( runlostlep ) {

		myContext.SetUseRl( true );
		myContext.SetJesDir( contextVars::kNominal );

		// Run lost lepton CR looper on all samples
		resetFile( crLostLep->GetPlotFileName() );
		resetFile( crLostLep->GetSystFileName() );
		for( sample* mySample : crLostLep->GetAllSamples() ) looperCR2lep( crLostLep, mySample );
	}

	////////////////////////////////////////////////
	// Run looperCR0b (the 0-bjet CR looper)

	if( run1lw ) {

		myContext.SetUseRl( false );
		myContext.SetJesDir( contextVars::kNominal );

		// Run 0-bjet CR looper on all samples
		resetFile( cr0bjets->GetPlotFileName() );
		resetFile( cr0bjets->GetSystFileName() );
		for( sample* mySample : cr0bjets->GetAllSamples() ) looperCR0b( cr0bjets, mySample );
	}

	/////////////////////////////////////////////////////
	// Run jet energy scale (JES) up and down variations

	if( runjes ) {

		// Reset output files
		resetFile( sr_jesup->GetPlotFileName());
		resetFile( cr2l_jesup->GetPlotFileName());
		resetFile( cr0b_jesup->GetPlotFileName());

		// Run ScanChain, looperCR2lep, and loopercr0b on JES up/down babies
		myContext.SetJesDir( contextVars::kUp );
		for( sample* mySample : sr_jesup->GetAllSamples()   ) ScanChain(    sr_jesup,   mySample );
		for( sample* mySample : cr2l_jesup->GetAllSamples() ) looperCR2lep( cr2l_jesup, mySample );
		for( sample* mySample : cr0b_jesup->GetAllSamples() ) looperCR0b(   cr0b_jesup, mySample );

		myContext.SetJesDir( contextVars::kDown );
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
		makeSignalEstimate( sr_signal, sig_genmet );
	}
	if( runcards )                makeDataCards( srAnalysis, sr_signal, crLostLep, cr0bjets );

	// Clean up /////////
	delete srAnalysis;
	delete crLostLep;
	delete cr0bjets;
	delete sr_signal;
	delete data;
	delete signal;
	delete tt2l;
	delete tt1l;
	delete singtop;
	delete wjets;
	delete dy;
	delete rare;

	return 0;
}


void resetFile( TString filename ) {
	TFile* f = new TFile( filename, "RECREATE" );
	delete f;
}
