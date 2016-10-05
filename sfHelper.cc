#include "sfHelper.h"
#include <iostream>


sfHelper myHelper; // Define an extern sfHelper



// Function definitions for "sfHelper" class

// Setup function
void sfHelper::Setup( bool is_fastsim, TH1D* counterHist, TH2F* nevtsHist=NULL, TH3D* counterHist_SMS=NULL )
{
	isFastsim = is_fastsim;
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

// Get a 7% upward/downward variation for the cr2l trigger efficiency
double sfHelper::Trig2lUp()   { return 1.07; }
double sfHelper::Trig2lDown() { return 0.93; }



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
	double Trig2lUp()      { return myHelper.Trig2lUp(); }
	double Trig2lDown()    { return myHelper.Trig2lDown(); }
}
