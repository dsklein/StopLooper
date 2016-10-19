#ifndef SFHELPER_H
#define SFHELPER_H


#include "TString.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3D.h"

#include "CMS3.h"


// SfHelper class:
// Provides a set of functions that will painlessly deliver you reweighting factors


class sfHelper {

public:

	sfHelper();
	~sfHelper();

	void Setup( bool is_fastsim, bool is_cr2l, TH1D* counterHist, TH2F* nevtsHist, TH3D* counterHist_SMS );

	double LepSFUp();
	double LepSFDown();
	double BtagHeavyUp();
	double BtagHeavyDown();
	double BtagLightUp();
	double BtagLightDown();
	double ISRUp();
	double ISRDown();
	double QSquaredUp();
	double QSquaredDown();
	double AlphaSUp();
	double AlphaSDown();
	double Unity();
	double TrigEff2l();
	double Trig2lUp();
	double Trig2lDown();
	double MetResSF();
	double MetResUp();
	double MetResDown();


private:
	bool isFastsim;
	bool isCR2l;
	TH1D* h_counter;
	TH2F* hist_nEvts;
	TH3D* h_counterSMS;
	TH2F* h_trigeff_cr2l;

	double lepnorm;
	double lepnorm_up;
	double lepnorm_down;
	double btagnorm;
	double btagnormHF_up;
	double btagnormLF_up;
	double btagnormHF_down;
	double btagnormLF_down;
	double qsquarednorm;
	double qsquarednorm_up;
	double qsquarednorm_down;
	double alphasnorm;
	double alphasnorm_up;
	double alphasnorm_down;

};


extern sfHelper myHelper;

namespace sfhelp {
	double LepSFUp();
	double LepSFDown();
	double BtagHeavyUp();
	double BtagHeavyDown();
	double BtagLightUp();
	double BtagLightDown();
	double ISRUp();
	double ISRDown();
	double QSquaredUp();
	double QSquaredDown();
	double AlphaSUp();
	double AlphaSDown();
	double Unity();
	double TrigEff2l();
	double Trig2lUp();
	double Trig2lDown();
	double MetResSF();
	double MetResUp();
	double MetResDown();
}


#endif
