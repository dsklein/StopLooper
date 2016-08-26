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

  void setup( bool is_fastsim, TH1D* counterHist, TH2F* nevtsHist, TH3D* counterHist_SMS );

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


private:
  bool isFastsim;
  TH1D* h_counter;
  TH2F* hist_nEvts;
  TH3D* h_counterSMS;

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
}


#endif
