// Class sfManager - a class to load scale factor histograms, and provide functions to calculate SFs for me


#include "sfManager.h"


// Constructor ///////////////////////////////

sfManager::sfManager( bool isFastsim, TString path, TH1D* counterHisto, variation var )
  : fastsim( isFastsim ),
	h_counter( counterHisto ),
	varType( var )
{

  bTagNormFactor = -999.99;
  TH1::SetDefaultSumw2();
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  // Load medium and veto electron SF histos
  TFile* f_el = new TFile( path + "/lepsf/kinematicBinSFele.root" );
  h_medEl_sf = (TH2D*) f_el->Get("CutBasedMedium")->Clone("medium_el_id");
  h_vetoEl_sf = (TH2D*) f_el->Get("CutBasedVeto")->Clone("veto_el_id");
  h_medEl_sf->Multiply( (TH2D*)f_el->Get("MiniIso0p1_vs_AbsEta") );
  h_vetoEl_sf->Multiply( (TH2D*)f_el->Get("MiniIso0p4_vs_AbsEta") );
  h_medEl_sf->SetDirectory(rootdir);
  h_vetoEl_sf->SetDirectory(rootdir);
  f_el->Close();

  // Load muon iso SF histo
  TFile* f_muIsotight = new TFile( path + "/lepsf/TnP_MuonID_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root" );
  TH2D* h_muIso_sf = (TH2D*) f_muIsotight->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_PF_pass_&_tag_IsoMu20_pass");

  // Load tight muon SF histo
  TFile* f_muIDtight = new TFile( path + "/lepsf/TnP_MuonID_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root" );
  h_tightMu_sf = (TH2D*) f_muIDtight->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_tag_IsoMu20_pass")->Clone("tight_mu_id");
  h_tightMu_sf->Multiply( h_muIso_sf );
  h_tightMu_sf->SetDirectory(rootdir);
  f_muIDtight->Close();

  // Load veto muon SF histo
  TFile* f_muIDveto = new TFile( path + "/lepsf/TnP_MuonID_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root" );
  h_vetoMu_sf = (TH2D*) f_muIDveto->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_tag_IsoMu20_pass")->Clone("veto_mu_id");
  h_vetoMu_sf->Multiply( h_muIso_sf );
  h_vetoMu_sf->SetDirectory(rootdir);
  f_muIDveto->Close();
  f_muIsotight->Close();

  // Load veto lepton efficiency histos
  TFile* f_vetoLepEff = new TFile( path + "/lepsf/lepeff__ttbar_powheg_pythia8_25ns__SRcuts.root" );
  h_vetoEl_eff = (TH2D*) f_vetoLepEff->Get("h2_lepEff_vetoSel_rebin_Eff_el")->Clone("veto_el_eff");
  h_vetoMu_eff = (TH2D*) f_vetoLepEff->Get("h2_lepEff_vetoSel_rebin_Eff_mu")->Clone("veto_mu_eff");
  h_vetoEl_eff->SetDirectory(rootdir);
  h_vetoMu_eff->SetDirectory(rootdir);
  f_vetoLepEff->Close();


  // Set up b-tag SF tool and readers
  calib = new BTagCalibration( "csvv2", (path + "/btagsf/data/run2_25ns/CSVv2.csv").Data() );
  calib_fastsim = new BTagCalibration( "CSV", (path + "/btagsf/data/run2_fastsim/CSV_13TEV_Combined_20_11_2015.csv").Data() );
  reader_heavy = new BTagCalibrationReader( calib, BTagEntry::OP_MEDIUM, "mujets", "central" );
  reader_light = new BTagCalibrationReader( calib, BTagEntry::OP_MEDIUM, "comb",   "central" );
  reader_fastsim = new BTagCalibrationReader( calib_fastsim, BTagEntry::OP_MEDIUM, "fastsim", "central" );

  // Get b-tagging efficiency histograms
  TFile* f_btagEff;
  if( fastsim ) f_btagEff = new TFile( path + "/btagsf/data/run2_fastsim/btageff__SMS-T1bbbb-T1qqqq_fastsim.root" );
  else          f_btagEff = new TFile( path + "/btagsf/data/run2_25ns/btageff__ttbar_powheg_pythia8_25ns.root" );
  h_btagEff_b = (TH2D*) f_btagEff->Get("h2_BTaggingEff_csv_med_Eff_b")->Clone("btag_eff_b");
  h_btagEff_c = (TH2D*) f_btagEff->Get("h2_BTaggingEff_csv_med_Eff_c")->Clone("btag_eff_c");
  h_btagEff_udsg = (TH2D*) f_btagEff->Get("h2_BTaggingEff_csv_med_Eff_udsg")->Clone("btag_eff_udsg");
  h_btagEff_b->SetDirectory(rootdir);
  h_btagEff_c->SetDirectory(rootdir);
  h_btagEff_udsg->SetDirectory(rootdir);
  f_btagEff->Close();

} // end constructor


// Destructor ////////////////////////////////

sfManager::~sfManager() {
  delete h_medEl_sf;
  delete h_vetoEl_sf;
  delete h_vetoEl_eff;
  delete h_tightMu_sf;
  delete h_vetoMu_sf;
  delete h_vetoMu_eff;

  delete calib;
  delete calib_fastsim;
  delete reader_heavy;
  delete reader_light;
  delete reader_fastsim;
}



//////////////////////////////////////////////
// Other functions ///////////////////////////

// Manually set the number of events in the sample
void sfManager::SetBtagNorm( double factor ) { bTagNormFactor = factor; }


// Reco electron SF
double sfManager::GetSF_el( double pt, double eta ) {
  double pt_bounded = std::max( 10., std::min( 99.9, pt ));
  double eta_bounded = std::min( 2.399, fabs(eta) );
  int bin = h_medEl_sf->FindBin( pt_bounded, eta_bounded );
  double sf = h_medEl_sf->GetBinContent( bin );
  if(      varType==kLepUp   ) return sf + h_medEl_sf->GetBinError( bin );
  else if( varType==kLepDown ) return sf - h_medEl_sf->GetBinError( bin );
  else return sf;
}

// Veto electron SF
double sfManager::GetSF_elVeto( double pt, double eta ) {
  double pt_bounded = std::max( 10., std::min( 99.9, pt ));
  double eta_bounded = std::min( 2.399, fabs(eta) );
  int bin = h_vetoEl_sf->FindBin( pt_bounded, eta_bounded );
  double sf = h_vetoEl_sf->GetBinContent( bin );
  if(      varType==kLepUp   ) return sf + h_vetoEl_sf->GetBinError( bin );
  else if( varType==kLepDown ) return sf - h_vetoEl_sf->GetBinError( bin );
  else return sf;
}

// Lost electron SF
double sfManager::GetSF_elLost( double pt, double eta ) {
  double pt_bounded = std::max( 10., std::min( 99.9, pt ));
  double eta_bounded = std::min( 2.399, fabs(eta) );
  double sf  = GetSF_elVeto( pt, eta );
  double eff = h_vetoEl_eff->GetBinContent( h_vetoEl_eff->FindBin( pt_bounded, eta_bounded ) );
  return ( 1.-sf*eff )/( 1.-eff);
}

// Reco muon SF
double sfManager::GetSF_mu( double pt, double eta ) {
  double pt_bounded = std::max( 10., std::min( 99.9, pt ));
  double eta_bounded = std::min( 2.399, fabs(eta) );
  int bin = h_tightMu_sf->FindBin( pt_bounded, eta_bounded );
  double sf = h_tightMu_sf->GetBinContent( bin );
  if(      varType==kLepUp   ) return sf + h_tightMu_sf->GetBinError( bin );
  else if( varType==kLepDown ) return sf - h_tightMu_sf->GetBinError( bin );
  else return sf;
}

// Veto muon SF
double sfManager::GetSF_muVeto( double pt, double eta ) {
  double pt_bounded = std::max( 10., std::min( 99.9, pt ));
  double eta_bounded = std::min( 2.399, fabs(eta) );
  int bin = h_vetoMu_sf->FindBin( pt_bounded, eta_bounded );
  double sf = h_vetoMu_sf->GetBinContent( bin );
  if(      varType==kLepUp   ) return sf + h_vetoMu_sf->GetBinError( bin );
  else if( varType==kLepDown ) return sf - h_vetoMu_sf->GetBinError( bin );
  else return sf;
}

// Lost muon SF
double sfManager::GetSF_muLost( double pt, double eta ) {
  double pt_bounded = std::max( 10., std::min( 99.9, pt ));
  double eta_bounded = std::min( 2.399, fabs(eta) );
  double sf  = GetSF_muVeto( pt_bounded, eta_bounded );
  double eff = h_vetoMu_eff->GetBinContent( h_vetoMu_eff->FindBin( pt_bounded, eta_bounded ) );
  return ( 1.-sf*eff )/( 1.-eff);
}

// Btagging SF (applied even to non-btagged and non-b-hadron jets)
double sfManager::GetSF_btag( double pt, double eta, int hadFlavor, double csv ) {

  BTagEntry::JetFlavor flavor;
  double sf = 1.;
  double eff = 1.;
  double weight = 1.;

  // Restrict our pt and eta to the ranges for which the btag SFs / efficiencies are valid
  double pt_sf  = std::max( 30.,   std::min(669.9,   pt) ); // 30 - 670
  double eta_sf = std::max( -2.399, std::min(2.399, eta) ); // -2.4 - 2.4
  double pt_eff  = std::max( 20.,   std::min(599.9,   pt) ); // 20 - 600
  double eta_eff = std::min( 2.799, fabs(eta) ); // 0 - 2.8
  double pt_fastsim = std::max( 20., std::min(799.9, pt) ); // 20 - 800

  // Choose which flavor to put in to the calibration reader
  if( hadFlavor==5 )        flavor = BTagEntry::FLAV_B;
  else if( hadFlavor == 4 ) flavor = BTagEntry::FLAV_C;
  else                      flavor = BTagEntry::FLAV_UDSG;

  // Get the bTag SF
  if( flavor == BTagEntry::FLAV_UDSG ) sf = reader_light->eval( flavor, eta_sf, pt_sf );
  else                                 sf = reader_heavy->eval( flavor, eta_sf, pt_sf );

  // If it's a fastsim sample, apply additional fastsim SF
  if( fastsim ) sf *= reader_fastsim->eval( flavor, eta_sf, pt_fastsim );

  // Get the b-tagging efficiency for the corresponding jet flavor
  if(      flavor == BTagEntry::FLAV_B ) eff = h_btagEff_b->GetBinContent( h_btagEff_b->FindBin(pt_eff, eta_eff) );
  else if( flavor == BTagEntry::FLAV_C ) eff = h_btagEff_c->GetBinContent( h_btagEff_c->FindBin(pt_eff, eta_eff) );
  else                                   eff = h_btagEff_udsg->GetBinContent( h_btagEff_udsg->FindBin(pt_eff, eta_eff) );

  // If it passes the b-tag threshold, apply b-tag SF.
  // If not, apply SF to the non-btag efficiency
  if( csv > 0.890 ) weight = sf;
  else              weight = ( 1. - eff*sf )/( 1. - eff );

  // Scale event weight to account for overall normalization of btag SFs
  // If this number has been provided manually, use it. Otherwise, derive it and store it.
  if( bTagNormFactor < 0. ) bTagNormFactor = (h_counter->GetBinContent(22)) / (h_counter->GetBinContent(14));
  weight *= bTagNormFactor;

  return weight;
}
