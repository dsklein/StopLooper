#include "sfHelper.h"
#include <iostream>


sfHelper myHelper; // Define an extern sfHelper



// Function definitions for "sfHelper" class

// Constructor and destructor
sfHelper::sfHelper()
{
	TFile trigeffFile( "reference-data/triggerefficiency_2lCR.root" );
	h_trigeff_cr2l = (TH2F*)trigeffFile.Get("twoDefficiencypass_gapsfilled")->Clone("trigeff_cr2l");
	h_trigeff_cr2l->SetDirectory(0);
	trigeffFile.Close();
}

sfHelper::~sfHelper()
{
	delete h_trigeff_cr2l;
}

// Setup function
void sfHelper::Setup( bool is_fastsim, bool is_cr2l, TH1D* counterHist, TH2F* nevtsHist=NULL, TH3D* counterHist_SMS=NULL )
{
	isFastsim = is_fastsim;
	isCR2l = is_cr2l;
	h_counter = counterHist;
	hist_nEvts = nevtsHist;
	h_counterSMS = counterHist_SMS;

	qsquarednorm    = h_counter->GetBinContent( 1 );
	alphasnorm      = h_counter->GetBinContent( 1 );
	qsquarednorm_up = h_counter->GetBinContent( 5 );
	qsquarednorm_down = h_counter->GetBinContent( 9 );
	alphasnorm_up   = h_counter->GetBinContent( 12 );
	alphasnorm_down = h_counter->GetBinContent( 13 );
	btagnorm        = h_counter->GetBinContent( 14 );
	btagnormHF_up   = h_counter->GetBinContent( 15 );
	btagnormLF_up   = h_counter->GetBinContent( 16 );
	btagnormHF_down = h_counter->GetBinContent( 17 );
	btagnormLF_down = h_counter->GetBinContent( 18 );
	lepnorm         = h_counter->GetBinContent( 28 ) * h_counter->GetBinContent( 31 ); // Lepton SF * veto lepton SF
	lepnorm_up      = h_counter->GetBinContent( 29 ) * h_counter->GetBinContent( 32 );
	lepnorm_down    = h_counter->GetBinContent( 30 ) * h_counter->GetBinContent( 33 );
}

// Get reweighting factor to vary the lepton and veto lepton SFs up
double sfHelper::LepSFUp() {
	double lepsf    = tas::weight_lepSF()    * tas::weight_vetoLepSF();
	double lepsf_up = tas::weight_lepSF_up() * tas::weight_vetoLepSF_up();
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		lepnorm    = h_counterSMS->GetBinContent( binx, biny, 27 ) * h_counterSMS->GetBinContent( binx, biny, 30 );
		lepnorm_up = h_counterSMS->GetBinContent( binx, biny, 28 ) * h_counterSMS->GetBinContent( binx, biny, 31);
	}
	return (lepsf_up / lepnorm_up) / (lepsf / lepnorm);
}

// Get reweighting factor to vary the lepton and veto lepton SFs down
double sfHelper::LepSFDown() {
	double lepsf      = tas::weight_lepSF()      * tas::weight_vetoLepSF();
	double lepsf_down = tas::weight_lepSF_down() * tas::weight_vetoLepSF_down();
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		lepnorm      = h_counterSMS->GetBinContent( binx, biny, 27 ) * h_counterSMS->GetBinContent( binx, biny, 30 );
		lepnorm_down = h_counterSMS->GetBinContent( binx, biny, 29 ) * h_counterSMS->GetBinContent( binx, biny, 32);;
	}
	return (lepsf_down / lepnorm_down) / (lepsf / lepnorm);
}

// Get reweighting factor to vary the btag heavy flavor SF up
double sfHelper::BtagHeavyUp() {
	double btagsf    = tas::weight_btagsf();
	double btagsf_up = tas::weight_btagsf_heavy_UP();
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		btagnorm      = h_counterSMS->GetBinContent( binx, biny, 14 );
		btagnormHF_up = h_counterSMS->GetBinContent( binx, biny, 15 );
	}
	return (btagsf_up / btagnormHF_up) / (btagsf / btagnorm);
}

// Get reweighting factor to vary the btag heavy flavor SF down
double sfHelper::BtagHeavyDown() {
	double btagsf      = tas::weight_btagsf();
	double btagsf_down = tas::weight_btagsf_heavy_DN();
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		btagnorm        = h_counterSMS->GetBinContent( binx, biny, 14 );
		btagnormHF_down = h_counterSMS->GetBinContent( binx, biny, 17 );
	}
	return (btagsf_down / btagnormHF_down) / (btagsf / btagnorm);
}

// Get reweighting factor to vary the btag light flavor SF up
double sfHelper::BtagLightUp() {
	double btagsf    = tas::weight_btagsf();
	double btagsf_up = tas::weight_btagsf_light_UP();
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		btagnorm      = h_counterSMS->GetBinContent( binx, biny, 14 );
		btagnormLF_up = h_counterSMS->GetBinContent( binx, biny, 16 );
	}
	return (btagsf_up / btagnormLF_up) / (btagsf / btagnorm);
}

// Get reweighting factor to vary the btag light flavor SF down
double sfHelper::BtagLightDown() {
	double btagsf      = tas::weight_btagsf();
	double btagsf_down = tas::weight_btagsf_light_DN();
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		btagnorm        = h_counterSMS->GetBinContent( binx, biny, 14 );
		btagnormLF_down = h_counterSMS->GetBinContent( binx, biny, 18 );
	}
	return (btagsf_down / btagnormLF_down) / (btagsf / btagnorm);
}

// Get reweighting factor to vary the ISR weight up
double sfHelper::ISRUp() {
	if( !isFastsim ) {
		std::cout << "Warning in sfHelper.cc: tried to apply ISR weight to a non-signal sample!" << std::endl;
		return 1.0;
	}
	double isr    = tas::weight_ISR();
	double isr_up = tas::weight_ISRup();
	int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
	int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
	double isrnorm    = h_counterSMS->GetBinContent( binx, biny, 19 );
	double isrnorm_up = h_counterSMS->GetBinContent( binx, biny, 20 );
	return (isr_up / isrnorm_up ) / (isr / isrnorm );
}

// Get reweighting factor to vary the ISR weight down
double sfHelper::ISRDown() {
	if( !isFastsim ) {
		std::cout << "Warning in sfHelper.cc: tried to apply ISR weight to a non-signal sample!" << std::endl;
		return 1.0;
	}
	double isr      = tas::weight_ISR();
	double isr_down = tas::weight_ISRdown();
	int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
	int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
	double isrnorm      = h_counterSMS->GetBinContent( binx, biny, 19 );
	double isrnorm_down = h_counterSMS->GetBinContent( binx, biny, 21 );
	return (isr_down / isrnorm_down ) / (isr / isrnorm );
}

// Get reweighting factor to vary Q^2 up
double sfHelper::QSquaredUp() {
	if( tas::genweights().size() < 111 ) return 1.;
	double qsquared    = tas::genweights().at(0);
	double qsquared_up = tas::genweights().at(4);
	if( qsquared < 0. || qsquared_up < 0. ) return 1.;
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		qsquarednorm    = h_counterSMS->GetBinContent( binx, biny, 1 );
		qsquarednorm_up = h_counterSMS->GetBinContent( binx, biny, 5 );
	}
	return (qsquared_up / qsquarednorm_up ) / (qsquared / qsquarednorm );
}

// Get reweighting factor to vary Q^2 down
double sfHelper::QSquaredDown() {
	if( tas::genweights().size() < 111 ) return 1.;
	double qsquared      = tas::genweights().at(0);
	double qsquared_down = tas::genweights().at(8);
	if( qsquared < 0. || qsquared_down < 0. ) return 1.;
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		qsquarednorm      = h_counterSMS->GetBinContent( binx, biny, 1 );
		qsquarednorm_down = h_counterSMS->GetBinContent( binx, biny, 9 );
	}
	return (qsquared_down / qsquarednorm_down ) / (qsquared / qsquarednorm );
}

// Get reweighting factor to vary alpha_s up
double sfHelper::AlphaSUp() {
	if( tas::genweights().size() < 111 ) return 1.;
	double alphas    = tas::genweights().at(0);
	double alphas_up = tas::genweights().at(109);
	if( alphas < 0. || alphas_up < 0. ) return 1.;
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		alphasnorm    = h_counterSMS->GetBinContent( binx, biny, 1 );
		alphasnorm_up = h_counterSMS->GetBinContent( binx, biny, 12 );
	}
	return (alphas_up / alphasnorm_up ) / (alphas / alphasnorm );
}

// Get reweighting factor to vary alpha_s down
double sfHelper::AlphaSDown() {
	if( tas::genweights().size() < 111 ) return 1.;
	double alphas      = tas::genweights().at(0);
	double alphas_down = tas::genweights().at(110);
	if( alphas < 0. || alphas_down < 0. ) return 1.;
	if( isFastsim ) {
		int binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
		int biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );
		alphasnorm      = h_counterSMS->GetBinContent( binx, biny, 1 );
		alphasnorm_down = h_counterSMS->GetBinContent( binx, biny, 13 );
	}
	return (alphas_down / alphasnorm_down ) / (alphas / alphasnorm );
}

// Dummy function that returns 1.0
double sfHelper::Unity() { return 1.0; }

// Get the 2-lepton CR trigger efficiency
double sfHelper::TrigEff2l() {
	double mymet = max(  250., min(499.99, double(tas::pfmet_rl()) ) );
	double myleppt = max( 20., min(499.99, double(tas::lep1_p4().pt()) ) );
	return h_trigeff_cr2l->GetBinContent( h_trigeff_cr2l->FindBin(mymet, myleppt) );
}

// Get reweighting factor to vary the cr2l trigger efficiency up by 7 percentage points
double sfHelper::Trig2lUp() {
	double oldeff = TrigEff2l();
	double neweff = oldeff + 0.07;
	return neweff / oldeff;
}

// Get reweighting factor to vary the cr2l trigger efficiency down by 7 percentage points
double sfHelper::Trig2lDown() {
	double oldeff = TrigEff2l();
	double neweff = oldeff - 0.07;
	return neweff / oldeff;
}

// Get MET resolution SF
double sfHelper::MetResSF() {
	double mymet = isCR2l ? tas::pfmet_rl() : tas::pfmet();
	double mymt2w = isCR2l ? tas::MT2W_rl() : tas::MT2W();

	if( tas::ngoodjets() == 2) {
		if(      mymet >= 450. ) return 0.776;
		else if( mymet >= 350. ) return 0.910;
		else if( mymet >= 250. ) return 1.048;
	}
	else if( tas::ngoodjets() == 3 ) {
		if(      mymet >= 550. ) return 0.658;
		else if( mymet >= 450. ) return 0.808;
		else if( mymet >= 350. ) return 0.928;
		else if( mymet >= 250. ) return 1.076;
	}
	else if( tas::ngoodjets() >= 4 ) {
		if( mymt2w < 200. && mymet >= 450. ) return 0.740;
		else if( mymet >= 650. ) return 0.564;
		else if( mymet >= 550. ) return 0.747;
		else if( mymet >= 450. ) return 0.873;
		else if( mymet >= 350. ) return 0.939;
		else if( mymet >= 250. ) return 1.076;
	}

	return 1.0;
}

// Get reweighting factor to vary the MET resolution SF up
double sfHelper::MetResUp() {
	double sf = MetResSF();
	double err = fabs(1.0-sf) / 2.;
	return (sf+err) / sf;
}

// Get reweighting factor to vary the MET resolution SF down
double sfHelper::MetResDown() {
	double sf = MetResSF();
	double err = fabs(1.0-sf) / 2.;
	return (sf-err) / sf;
}



namespace sfhelp {
	double LepSFUp()       { return myHelper.LepSFUp(); }
	double LepSFDown()     { return myHelper.LepSFDown(); }
	double BtagHeavyUp()   { return myHelper.BtagHeavyUp(); }
	double BtagHeavyDown() { return myHelper.BtagHeavyDown(); }
	double BtagLightUp()   { return myHelper.BtagLightUp(); }
	double BtagLightDown() { return myHelper.BtagLightDown(); }
	double ISRUp()         { return myHelper.ISRUp(); }
	double ISRDown()       { return myHelper.ISRDown(); }
	double QSquaredUp()    { return myHelper.QSquaredUp(); }
	double QSquaredDown()  { return myHelper.QSquaredDown(); }
	double AlphaSUp()      { return myHelper.AlphaSUp(); }
	double AlphaSDown()    { return myHelper.AlphaSDown(); }
	double Unity()         { return myHelper.Unity(); }
	double TrigEff2l()     { return myHelper.TrigEff2l(); }
	double Trig2lUp()      { return myHelper.Trig2lUp(); }
	double Trig2lDown()    { return myHelper.Trig2lDown(); }
	double MetResSF()      { return myHelper.MetResSF(); }
	double MetResUp()      { return myHelper.MetResUp(); }
	double MetResDown()    { return myHelper.MetResDown(); }
}
