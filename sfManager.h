// Class sfManager - a class to load scale factor histograms, and provide functions to calculate SFs for me

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"

#include "btagsf/BTagCalibrationStandalone.h"


class sfManager{

public:
  sfManager( bool isFastsim, TString path, TH1D* counterHisto );
  ~sfManager();

  double GetSF_el( double pt, double eta );
  double GetSF_elVeto( double pt, double eta );
  double GetSF_elLost( double pt, double eta );
  double GetSF_mu( double pt, double eta );
  double GetSF_muVeto( double pt, double eta );
  double GetSF_muLost( double pt, double eta );

  double GetSF_btag( double pt, double eta, int hadFlavor, double csv );

  void SetBtagNorm( double factor );

private:
  bool fastsim;
  double bTagNormFactor;

  TH1D* h_counter;
  TH2D* h_medEl_sf;
  TH2D* h_vetoEl_sf;
  TH2D* h_vetoEl_eff;
  TH2D* h_tightMu_sf;
  TH2D* h_vetoMu_sf;
  TH2D* h_vetoMu_eff;
  TH2D* h_btagEff_b;
  TH2D* h_btagEff_c;
  TH2D* h_btagEff_udsg;

  BTagCalibration* calib;
  BTagCalibration* calib_fastsim;
  BTagCalibrationReader* reader_heavy;
  BTagCalibrationReader* reader_light;
  BTagCalibrationReader* reader_fastsim;

};
