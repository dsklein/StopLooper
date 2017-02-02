#ifndef SFHELPER_H
#define SFHELPER_H

#include "TEfficiency.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "Math/VectorUtil.h"

#include "CMS3.h"

#include "contextVars.h"


// SfHelper class:
// Provides a set of functions that will painlessly deliver you reweighting factors


class sfHelper {

public:

	sfHelper();
	~sfHelper();

	void Setup( bool is_fastsim, TH1D* counterHist, TH2F* nevtsHist, TH3D* counterHist_SMS );
	void PrepSignal();
	void SetCorridor( bool corridor );

	double LepSF();
	double LepSFUp();
	double LepSFDown();
	double LepSFfastsim();
	double LepSFfastsimUp();
	double LepSFfastsimDown();
	double BtagSF();
	double BtagHeavyUp();
	double BtagHeavyDown();
	double BtagLightUp();
	double BtagLightDown();
	double BtagFSUp();
	double BtagFSDown();
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
	double MetResSF_corr();
	double MetResCorrectionCorridor();
	double MetResUp();
	double MetResDown();
	double TopSystPtSF();
	double TopSystPtUp();
	double TopSystPtDown();
	double Contam1lwUp();
	double Contam1lwDown();
	double ISRnJetsSF();
	double ISRnJetsUp();
	double ISRnJetsDown();
	double LumiUp();
	double PDFUp();
	double PDFDown();
	double StopXsecUp();
	double StopXsecDown();

private:
	bool isFastsim;
	int topptBin;
	int nEvts;
	int binx;
	int biny;
	TH1D* h_counter;
	TH2F* hist_nEvts;
	TH3D* h_counterSMS;
	TEfficiency* eff_1l_el;
	TEfficiency* eff_1l_mu;
	TEfficiency* eff_2l_el;
	TEfficiency* eff_2l_mu;
	double eff_err_up;
	double eff_err_down;
	TH1D* h_sf_toppt;
	bool isCorridor;
	TH1D* h_stop_xsec;

	double lepnorm;
	double lepnorm_up;
	double lepnorm_down;
	double lepnorm_fastsim;
	double lepnorm_fs_up;
	double lepnorm_fs_down;
	double btagnorm;
	double btagnormHF_up;
	double btagnormLF_up;
	double btagnormHF_down;
	double btagnormLF_down;
	double isrnorm;
	double isrnorm_up;
	double isrnorm_down;
	double btagnormFS_up;
	double btagnormFS_down;
	double isrnjetsnorm;
	double isrnjetsnorm_up;
	double isrnjetsnorm_down;
	double qsquarednorm;
	double qsquarednorm_up;
	double qsquarednorm_down;
	double PDFnorm_up;
	double PDFnorm_down;
	double alphasnorm;
	double alphasnorm_up;
	double alphasnorm_down;
};


extern sfHelper myHelper;

namespace sfhelp {
	double LepSF();
	double LepSFUp();
	double LepSFDown();
	double LepSFfastsim();
	double LepSFfastsimUp();
	double LepSFfastsimDown();
	double BtagSF();
	double BtagHeavyUp();
	double BtagHeavyDown();
	double BtagLightUp();
	double BtagLightDown();
	double BtagFSUp();
	double BtagFSDown();
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
	double MetResSF_corr();
	double MetResCorrectionCorridor();
	double MetResUp();
	double MetResDown();
	double TopSystPtSF();
	double TopSystPtUp();
	double TopSystPtDown();
	double Contam1lwUp();
	double Contam1lwDown();
	double ISRnJetsSF();
	double ISRnJetsUp();
	double ISRnJetsDown();
	double LumiUp();
	double PDFUp();
	double PDFDown();
	double StopXsecUp();
	double StopXsecDown();
}


#endif
