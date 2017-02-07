#include "sfHelper.h"
#include <iostream>


sfHelper myHelper; // Define an extern sfHelper

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;


// Function definitions for "sfHelper" class

// Constructor and destructor
sfHelper::sfHelper()
{
	TFile trigeffFile( "reference-files/TriggerEff.root", "READ" );
	eff_1l_el = (TEfficiency*)trigeffFile.Get("Efficiency_ge1l_el")->Clone("trigeff_1l_el");
	eff_1l_mu = (TEfficiency*)trigeffFile.Get("Efficiency_ge1l_mu")->Clone("trigeff_1l_mu");
	eff_2l_el = (TEfficiency*)trigeffFile.Get("Efficiency_ge2l_metrl_el")->Clone("trigeff_2l_el");
	eff_2l_mu = (TEfficiency*)trigeffFile.Get("Efficiency_ge2l_metrl_mu")->Clone("trigeff_2l_mu");
	eff_1l_el->SetDirectory(0);
	eff_1l_mu->SetDirectory(0);
	eff_2l_el->SetDirectory(0);
	eff_2l_mu->SetDirectory(0);
	trigeffFile.Close();

	TFile topptFile( "reference-files/sf_top_system_pt.root", "READ" );
	h_sf_toppt = (TH1D*)topptFile.Get("topsyst_pt_sf")->Clone("top_system_pt_sf");
	h_sf_toppt->SetDirectory(0);
	topptFile.Close();

	TFile stopxsecFile( "reference-files/xsec_stop_13TeV.root", "READ" );
	h_stop_xsec = (TH1D*)stopxsecFile.Get("stop")->Clone("stop_xsec");
	h_stop_xsec->SetDirectory(0);
	stopxsecFile.Close();
}

sfHelper::~sfHelper()
{
	delete eff_1l_el;
	delete eff_1l_mu;
	delete eff_2l_el;
	delete eff_2l_mu;
	delete h_sf_toppt;
	delete h_stop_xsec;
}

// Setup function
void sfHelper::Setup( bool is_fastsim, TH1D* counterHist=NULL, TH2F* nevtsHist=NULL, TH3D* counterHist_SMS=NULL )
{
	isFastsim = is_fastsim;
	isCorridor = false;
	h_counter = counterHist;
	hist_nEvts = nevtsHist;
	h_counterSMS = counterHist_SMS;
	topptBin = -99;

	TString filename = TString( gFile->GetName() );
	if(      filename.Contains("data_") ) sampleType = kData;
	else if( filename.Contains("Signal") ) sampleType = kSignal;
	else if( filename.Contains("ttbar_diLept") ) sampleType = ktt2l;
	else if( filename.Contains("ttbar_singleLept") ) sampleType = ktt1l;
	else if( filename.Contains("ttWJets") ) sampleType = kttW;
	else if( filename.Contains("ttZJets") ) sampleType = kttZ;
	else if( filename.Contains("ttbar") ) sampleType = kttbar;
	else if( filename.Contains("WWTo") ||
	         filename.Contains("WZTo") ||
	         filename.Contains("ZZTo") ) sampleType = kDiboson;
	else if( filename.Contains("W_5f_powheg_pythia8") ) sampleType = ktW;
	else if( filename.Contains("ch_4f_") ) sampleType = kSingletop;
	else if( filename.Contains("JetsToLNu_") ) sampleType = kWjets;
	else sampleType = kOther;

	if( !h_counter ) return;

	nEvts = h_counter->GetBinContent( 22 );

	qsquarednorm    = h_counter->GetBinContent( 1 );
	alphasnorm      = h_counter->GetBinContent( 1 );
	qsquarednorm_up = h_counter->GetBinContent( 5 );
	qsquarednorm_down = h_counter->GetBinContent( 9 );
	PDFnorm_up      = h_counter->GetBinContent( 10 );
	PDFnorm_down    = h_counter->GetBinContent( 11 );
	alphasnorm_up   = h_counter->GetBinContent( 12 );
	alphasnorm_down = h_counter->GetBinContent( 13 );
	btagnorm        = h_counter->GetBinContent( 14 );
	btagnormHF_up   = h_counter->GetBinContent( 15 );
	btagnormLF_up   = h_counter->GetBinContent( 16 );
	btagnormHF_down = h_counter->GetBinContent( 17 );
	btagnormLF_down = h_counter->GetBinContent( 18 );
	btagnormFS_up   = h_counter->GetBinContent( 23 );
	btagnormFS_down = h_counter->GetBinContent( 24 );
	isrnjetsnorm    = h_counter->GetBinContent( 25 );
	isrnjetsnorm_up = h_counter->GetBinContent( 26 );
	isrnjetsnorm_down = h_counter->GetBinContent( 27 );
	lepnorm         = h_counter->GetBinContent( 28 ) * h_counter->GetBinContent( 31 ); // Lepton SF * veto lepton SF
	lepnorm_up      = h_counter->GetBinContent( 29 ) * h_counter->GetBinContent( 32 );
	lepnorm_down    = h_counter->GetBinContent( 30 ) * h_counter->GetBinContent( 33 );
	lepnorm_fastsim = h_counter->GetBinContent( 34 );
	lepnorm_fs_up   = h_counter->GetBinContent( 35 );
	lepnorm_fs_down = h_counter->GetBinContent( 36 );
}

// Special setup function for SUSY signal samples
void sfHelper::PrepSignal() {
	if( !hist_nEvts ) return;
	if( !h_counterSMS ) return;

	isCorridor = false;

	binx = h_counterSMS->GetXaxis()->FindBin( tas::mass_stop() );
	biny = h_counterSMS->GetYaxis()->FindBin( tas::mass_lsp() );

	nEvts = hist_nEvts->GetBinContent( binx, biny );

	qsquarednorm    = h_counterSMS->GetBinContent( binx, biny, 1 );
	alphasnorm      = h_counterSMS->GetBinContent( binx, biny, 1 );
	qsquarednorm_up = h_counterSMS->GetBinContent( binx, biny, 5 );
	qsquarednorm_down = h_counterSMS->GetBinContent( binx, biny, 9 );
	PDFnorm_up      = h_counterSMS->GetBinContent( binx, biny, 10 );
	PDFnorm_down    = h_counterSMS->GetBinContent( binx, biny, 11 );
	alphasnorm_up   = h_counterSMS->GetBinContent( binx, biny, 12 );
	alphasnorm_down = h_counterSMS->GetBinContent( binx, biny, 13 );
	btagnorm        = h_counterSMS->GetBinContent( binx, biny, 14 );
	btagnormHF_up   = h_counterSMS->GetBinContent( binx, biny, 15 );
	btagnormLF_up   = h_counterSMS->GetBinContent( binx, biny, 16 );
	btagnormHF_down = h_counterSMS->GetBinContent( binx, biny, 17 );
	btagnormLF_down = h_counterSMS->GetBinContent( binx, biny, 18 );
	isrnorm         = h_counterSMS->GetBinContent( binx, biny, 19 );
	isrnorm_up      = h_counterSMS->GetBinContent( binx, biny, 20 );
	isrnorm_down    = h_counterSMS->GetBinContent( binx, biny, 21 );
	btagnormFS_up   = h_counterSMS->GetBinContent( binx, biny, 22 );
	btagnormFS_down = h_counterSMS->GetBinContent( binx, biny, 23 );
	isrnjetsnorm    = h_counterSMS->GetBinContent( binx, biny, 24 );
	isrnjetsnorm_up = h_counterSMS->GetBinContent( binx, biny, 25 );
	isrnjetsnorm_down = h_counterSMS->GetBinContent( binx, biny, 26 );
	lepnorm         = h_counterSMS->GetBinContent( binx, biny, 27 ) * h_counterSMS->GetBinContent( binx, biny, 30 ); // Lepton SF * veto lepton SF
	lepnorm_up      = h_counterSMS->GetBinContent( binx, biny, 28 ) * h_counterSMS->GetBinContent( binx, biny, 31 );
	lepnorm_down    = h_counterSMS->GetBinContent( binx, biny, 29 ) * h_counterSMS->GetBinContent( binx, biny, 32 );
	lepnorm_fastsim = h_counterSMS->GetBinContent( binx, biny, 33 );
	lepnorm_fs_up   = h_counterSMS->GetBinContent( binx, biny, 34 );
	lepnorm_fs_down = h_counterSMS->GetBinContent( binx, biny, 35 );
}

void sfHelper::SetCorridor( bool corridor ) {	isCorridor = corridor; }

// Get weight for lepton SF
double sfHelper::LepSF() {
	double sf = tas::weight_lepSF() * tas::weight_vetoLepSF();
	return sf * nEvts * nEvts / lepnorm; // Two factors of nEvts for lepton and vetolepton normalization
}

// Get reweighting factor to vary the lepton and veto lepton SFs up
double sfHelper::LepSFUp() {
	double lepsf    = tas::weight_lepSF()    * tas::weight_vetoLepSF();
	double lepsf_up = tas::weight_lepSF_up() * tas::weight_vetoLepSF_up();
	return (lepsf_up / lepnorm_up) / (lepsf / lepnorm);
}

// Get reweighting factor to vary the lepton and veto lepton SFs down
double sfHelper::LepSFDown() {
	double lepsf      = tas::weight_lepSF()      * tas::weight_vetoLepSF();
	double lepsf_down = tas::weight_lepSF_down() * tas::weight_vetoLepSF_down();
	return (lepsf_down / lepnorm_down) / (lepsf / lepnorm);
}

// Get weight for fastsim lepton SF
double sfHelper::LepSFfastsim() {
	if( !isFastsim ) return 1.0;
	double sf = tas::weight_lepSF_fastSim();
	return sf * nEvts / lepnorm_fastsim;
}

// Get reweighting factor to vary the fastsim lepton SF up
double sfHelper::LepSFfastsimUp() {
	if( !isFastsim ) return 1.0;
	double sf    = tas::weight_lepSF_fastSim();
	double sf_up = tas::weight_lepSF_fastSim_up();
	return (sf_up / lepnorm_fs_up) / (sf / lepnorm_fastsim);
}

// Get reweighting factor to vary the fastsim lepton SF down
double sfHelper::LepSFfastsimDown() {
	if( !isFastsim ) return 1.0;
	double sf      = tas::weight_lepSF_fastSim();
	double sf_down = tas::weight_lepSF_fastSim_down();
	return (sf_down / lepnorm_fs_down) / (sf / lepnorm_fastsim);
}

// Get weight for b-tag SF
double sfHelper::BtagSF() {
	double sf = tas::weight_btagsf();
	return sf * nEvts / btagnorm;
}

// Get reweighting factor to vary the btag heavy flavor SF up
double sfHelper::BtagHeavyUp() {
	double btagsf    = tas::weight_btagsf();
	double btagsf_up = tas::weight_btagsf_heavy_UP();
	return (btagsf_up / btagnormHF_up) / (btagsf / btagnorm);
}

// Get reweighting factor to vary the btag heavy flavor SF down
double sfHelper::BtagHeavyDown() {
	double btagsf      = tas::weight_btagsf();
	double btagsf_down = tas::weight_btagsf_heavy_DN();
	return (btagsf_down / btagnormHF_down) / (btagsf / btagnorm);
}

// Get reweighting factor to vary the btag light flavor SF up
double sfHelper::BtagLightUp() {
	double btagsf    = tas::weight_btagsf();
	double btagsf_up = tas::weight_btagsf_light_UP();
	return (btagsf_up / btagnormLF_up) / (btagsf / btagnorm);
}

// Get reweighting factor to vary the btag light flavor SF down
double sfHelper::BtagLightDown() {
	double btagsf      = tas::weight_btagsf();
	double btagsf_down = tas::weight_btagsf_light_DN();
	return (btagsf_down / btagnormLF_down) / (btagsf / btagnorm);
}

// Get reweighting factor to vary the btag fastsim SF up
double sfHelper::BtagFSUp() {
	if( !isFastsim ) return 1.0;
	double sf = tas::weight_btagsf_fastsim_UP();
	return sf * nEvts / btagnormFS_up;
}

// Get reweighting factor to vary the btag fastsim SF down
double sfHelper::BtagFSDown() {
	if( !isFastsim ) return 1.0;
	double sf = tas::weight_btagsf_fastsim_DN();
	return sf * nEvts / btagnormFS_down;
}

// Get reweighting factor to vary the ISR weight up
double sfHelper::ISRUp() {
	if( !isFastsim ) return 1.0;
	double isr    = tas::weight_ISR();
	double isr_up = tas::weight_ISRup();
	return (isr_up / isrnorm_up ) / (isr / isrnorm );
}

// Get reweighting factor to vary the ISR weight down
double sfHelper::ISRDown() {
	if( !isFastsim ) return 1.0;
	double isr      = tas::weight_ISR();
	double isr_down = tas::weight_ISRdown();
	return (isr_down / isrnorm_down ) / (isr / isrnorm );
}

// Get reweighting factor to vary Q^2 up
double sfHelper::QSquaredUp() {
	if( tas::genweights().size() < 110 ) return 1.;
	double qsquared    = tas::genweights().at(0);
	double qsquared_up = tas::genweights().at(4);
	if( qsquared < 0. || qsquared_up < 0. ) return 1.;
	return (qsquared_up / qsquarednorm_up ) / (qsquared / qsquarednorm );
}

// Get reweighting factor to vary Q^2 down
double sfHelper::QSquaredDown() {
	if( tas::genweights().size() < 110 ) return 1.;
	double qsquared      = tas::genweights().at(0);
	double qsquared_down = tas::genweights().at(8);
	if( qsquared < 0. || qsquared_down < 0. ) return 1.;
	return (qsquared_down / qsquarednorm_down ) / (qsquared / qsquarednorm );
}

// Get reweighting factor to vary alpha_s up
double sfHelper::AlphaSUp() {
	if( tas::genweights().size() < 110 ) return 1.;
	double alphas    = tas::genweights().at(0);
	double alphas_up = tas::genweights().at(109);
	if( alphas < 0. || alphas_up < 0. ) return 1.;
	return (alphas_up / alphasnorm_up ) / (alphas / alphasnorm );
}

// Get reweighting factor to vary alpha_s down
double sfHelper::AlphaSDown() {
	if( tas::genweights().size() < 110 ) return 1.;
	double alphas      = tas::genweights().at(0);
	double alphas_down = tas::genweights().at(110);
	if( alphas < 0. || alphas_down < 0. ) return 1.;
	return (alphas_down / alphasnorm_down ) / (alphas / alphasnorm );
}

// Dummy function that returns 1.0
double sfHelper::Unity() { return 1.0; }

// Get the 2-lepton CR trigger efficiency
double sfHelper::TrigEff2l() {
	double mymet = max(  150., min(499.99, double(context::Met()) ) );
	double myleppt = max( 15., min( 49.99, double(tas::lep1_p4().pt()) ) );

	TEfficiency* thisEff;
	if(      myContext.GetUseRl() && abs(tas::lep1_pdgid())==11 ) thisEff = eff_2l_el;
	else if( myContext.GetUseRl() && abs(tas::lep1_pdgid())==13 ) thisEff = eff_2l_mu;
	else if( abs(tas::lep1_pdgid())==11 ) thisEff = eff_1l_el;
	else if( abs(tas::lep1_pdgid())==13 ) thisEff = eff_1l_mu;

	int binnum = thisEff->FindFixBin( myleppt, mymet );
	double eff = thisEff->GetEfficiency( binnum );
	eff_err_up = thisEff->GetEfficiencyErrorUp( binnum );
	eff_err_down = thisEff->GetEfficiencyErrorLow( binnum );

	if( !myContext.GetUseRl() ) return 1.0; // HJ recommends we only use errors, and not SF, in 1-lepton regions
	return eff;
}

// Get reweighting factor to vary the cr2l trigger efficiency up by 7 percentage points
double sfHelper::Trig2lUp() {
	double oldeff = TrigEff2l();
	double neweff = oldeff + eff_err_up;
	return neweff / oldeff;
}

// Get reweighting factor to vary the cr2l trigger efficiency down by 7 percentage points
double sfHelper::Trig2lDown() {
	double oldeff = TrigEff2l();
	double neweff = oldeff - eff_err_down;
	return neweff / oldeff;
}

// Get MET resolution SF
double sfHelper::MetResSF() {
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar &&
	    sampleType != kWjets &&
	    sampleType != ktW ) return 1.0;

	int njets = context::ngoodjets();
	double met = context::Met();
	double modtop = context::TopnessMod();
	int ntightbs = context::ntightbtags();
	double mlb = context::Mlb_closestb();
	if( ntightbs == 0 && // If we're in CR0b, use leading CSV jet
	    (context::Mlb_closestb() >= 175. || context::ngoodbtags() == 0 ) ) mlb = context::Mlb_lead_bdiscr();

	if( njets < 4 && modtop >= 10. ) {
		if( mlb < 175. ) {
			if(      met >= 600. ) return 1.24; // Region A
			else if( met >= 450. ) return 1.14;
			else if( met >= 350. ) return 1.07;
			else if( met >= 250. ) return 0.98;
		}
		else if( mlb >= 175. && ntightbs >= 1 ) {
			if(      met >= 600. ) return 1.24; // B
			else if( met >= 450. ) return 1.14;
			else if( met >= 250. ) return 0.99;
		}
	}
	else if( njets >= 4 ) {
		if( modtop < 0. ) {
			if( mlb < 175. ) {
				if(      met >= 650. ) return 1.18; // C
				else if( met >= 550. ) return 1.21;
				else if( met >= 450. ) return 1.12;
				else if( met >= 350. ) return 1.06;
				else if( met >= 250. ) return 0.98;
			}
			else if( mlb >= 175. && ntightbs >= 1 ) {
				if(      met >= 550. ) return 1.20; // D
				else if( met >= 450. ) return 1.12;
				else if( met >= 350. ) return 1.06;
				else if( met >= 250. ) return 0.98;
			}
		}
		else if( modtop < 10. ) {
			if( mlb < 175. ) {
				if(      met >= 550. ) return 1.13; // E
				else if( met >= 350. ) return 1.08;
				else if( met >= 250. ) return 0.97;
			}
			else if( mlb >= 175. && ntightbs >= 1 ) {
				if(      met >= 450. ) return 1.12; // F
				else if( met >= 250. ) return 0.99;
			}
		}
		else if( modtop >= 10. ) {
			if( mlb < 175. ) {
				if(      met >= 600. ) return 1.11; // F
				else if( met >= 450. ) return 1.12;
				else if( met >= 350. ) return 1.08;
				else if( met >= 250. ) return 0.97;
			}
			else if( mlb >= 175. && ntightbs >= 1 ) {
				if(      met >= 450. ) return 1.12; // H
				else if( met >= 250. ) return 0.99;
			}
		}
	}

	return 1.0;
}

// MET resolution SF for corridor regions
double sfHelper::MetResSF_corr() {
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar &&
	    sampleType != kWjets &&
	    sampleType != ktW ) return 1.0;

	double met = context::Met();
	if(      met >= 550. ) return 1.14;
	else if( met >= 450. ) return 1.14;
	else if( met >= 350. ) return 1.05;
	else if( met >= 250. ) return 0.97;

	return 1.0;
}

// Get the correction factor to convert from Mlb-binned to corridor SF values
double sfHelper::MetResCorrectionCorridor() { return MetResSF_corr() / MetResSF(); }

// Get reweighting factor to vary the MET resolution SF up
double sfHelper::MetResUp() {
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar &&
	    sampleType != kWjets &&
	    sampleType != ktW ) return 1.0;

	double sf = isCorridor ? MetResSF_corr() : MetResSF();
	double err = fabs(1.0-sf) / 2.;
	return (sf+err) / sf;
}

// Get reweighting factor to vary the MET resolution SF down
double sfHelper::MetResDown() {
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar &&
	    sampleType != kWjets &&
	    sampleType != ktW ) return 1.0;

	double sf = isCorridor ? MetResSF_corr() : MetResSF();
	double err = fabs(1.0-sf) / 2.;
	return (sf+err) / sf;
}

// Get Top system pT scale factor
double sfHelper::TopSystPtSF() {
	topptBin = -99;
	if( !tas::is2lep() ) return 1.0;
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar &&
	    sampleType != ktW ) return 1.0;

	LorentzVector topsystem_p4( 0., 0., 0., 0. );

	topsystem_p4 += tas::lep1_p4();
	if( tas::nvetoleps() > 1 ) topsystem_p4 += tas::lep2_p4();

	int idx1 = -1;
	int idx2 = -1;
	double csv1 = -99.9;
	double csv2 = -99.9;
	double csv_new;

	// Find the two highest-CSV jet indices
	for( uint i=0; i<context::ak4pfjets_CSV().size(); i++ ) {
		csv_new = context::ak4pfjets_CSV().at(i);
		if( csv_new > csv1 ) {
			csv2 = csv1;
			idx2 = idx1;
			csv1 = csv_new;
			idx1 = i;
		}
		else if( csv_new > csv2 ) {
			csv2 = csv_new;
			idx2 = i;
		}
	}
	if( idx1 >= 0 ) topsystem_p4 += context::ak4pfjets_p4().at(idx1);
	if( idx2 >= 0 ) topsystem_p4 += context::ak4pfjets_p4().at(idx2);

	// Use MET without 2nd lep added
	double myMet = context::Met_no2ndlep();
	double myMetPhi = context::MetPhi_no2ndlep();
	LorentzVector met_p4( myMet*cos(myMetPhi), myMet*sin(myMetPhi), 0.0, myMet );
	topsystem_p4 += met_p4;

	double system_pt = topsystem_p4.Pt();
	topptBin = h_sf_toppt->FindBin(system_pt);
	return h_sf_toppt->GetBinContent(topptBin);
}

// Get reweighting factor to vary the top system pT SF up
double sfHelper::TopSystPtUp() {
	if( !tas::is2lep() ) return 1.0;
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar ) return 1.0;
	if( topptBin < 0 ) return 1.0;

	double sf  = h_sf_toppt->GetBinContent(topptBin);
	double err = h_sf_toppt->GetBinError(topptBin);
	return (sf+err) / sf;
}

// Get reweighting factor to vary the top system pT SF down
double sfHelper::TopSystPtDown() {
	if( !tas::is2lep() ) return 1.0;
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar ) return 1.0;
	if( topptBin < 0 ) return 1.0;

	double sf  = h_sf_toppt->GetBinContent(topptBin);
	double err = h_sf_toppt->GetBinError(topptBin);
	return (sf-err) / sf;
}

// Get reiweighting factor to vary the non-1l-from-W CR contamination up
double sfHelper::Contam1lwUp() {
	if(      tas::isZtoNuNu()     ) return 1.5;
	else if( tas::is2lep()        ) return 1.5;
	else if( tas::is1lepFromTop() ) return 1.5;
	return 1.0;
}

// Get reiweighting factor to vary the non-1l-from-W CR contamination down
double sfHelper::Contam1lwDown() {
	if(      tas::isZtoNuNu()     ) return 0.5;
	else if( tas::is2lep()        ) return 0.5;
	else if( tas::is1lepFromTop() ) return 0.5;
	return 1.0;
}

// Get ISR NJets scale factor
double sfHelper::ISRnJetsSF() {
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar &&
	    sampleType != kSignal ) return 1.0;

	double isrnjsf = tas::weight_ISRnjets();
	return isrnjsf * (nEvts / isrnjetsnorm);
}

// Get reweighting factor to vary the ISR NJets scale factor up
double sfHelper::ISRnJetsUp() {
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar &&
	    sampleType != kSignal ) return 1.0;

	double isrnjsf    = tas::weight_ISRnjets();
	double isrnjsf_up = tas::weight_ISRnjets_UP();
	return (isrnjsf_up / isrnjetsnorm_up) / (isrnjsf / isrnjetsnorm );
}

// Get reweighting factor to vary the ISR NJets scale factor down
double sfHelper::ISRnJetsDown() {
	if( sampleType != ktt2l &&
	    sampleType != ktt1l &&
	    sampleType != kttbar &&
	    sampleType != kSignal ) return 1.0;

	double isrnjsf      = tas::weight_ISRnjets();
	double isrnjsf_down = tas::weight_ISRnjets_DN();
	return (isrnjsf_down / isrnjetsnorm_down) / (isrnjsf / isrnjetsnorm );
}

// Get flat 6.2% uncertainty on lumi
double sfHelper::LumiUp() {
	if( tas::is_data() ) return 1.0;
	return 1.062;
}

// Get reweighting factor to vary PDF up
double sfHelper::PDFUp() {
	if( tas::genweights().size() < 110 ) return 1.;
	if( PDFnorm_up <= 0. ) return 1.;
	return fabs( tas::pdf_up_weight() * nEvts / PDFnorm_up );
}

// Get reweighting factor to vary PDF down
double sfHelper::PDFDown() {
	if( tas::genweights().size() < 110 ) return 1.;
	if( PDFnorm_down <= 0. ) return 1.;
	return fabs( tas::pdf_down_weight() * nEvts / PDFnorm_down );
}

// Get reweighting factor to vary W+HF xsec up
double sfHelper::WhfXsecUp() {
	if( sampleType != kWjets ) return 1.0;
	bool hasHFjet = false;

	for( int flavor : tas::ak4pfjets_hadron_flavor() ) {
		if( abs(flavor) == 5 ||
		    abs(flavor) == 4 ) {
			hasHFjet = true;
			break;
		}
	}

	if( hasHFjet ) return 1.5;
	return 1.0;
}

// Get reweighting factor to vary W+HF xsec down
double sfHelper::WhfXsecDown() {
	if( sampleType != kWjets ) return 1.0;
	bool hasHFjet = false;

	for( int flavor : tas::ak4pfjets_hadron_flavor() ) {
		if( abs(flavor) == 5 ||
		    abs(flavor) == 4 ) {
			hasHFjet = true;
			break;
		}
	}

	if( hasHFjet ) return 0.5;
	return 1.0;
}

// Get rewieighting factor to vary the stop xsec up
double sfHelper::StopXsecUp() {
	int bin = h_stop_xsec->FindBin( tas::mass_stop() );
	double xsec = h_stop_xsec->GetBinContent( bin );
	double err  = h_stop_xsec->GetBinError(   bin );
	return xsec + err;
}

// Get rewieighting factor to vary the stop xsec down
double sfHelper::StopXsecDown() {
	int bin = h_stop_xsec->FindBin( tas::mass_stop() );
	double xsec = h_stop_xsec->GetBinContent( bin );
	double err  = h_stop_xsec->GetBinError(   bin );
	return xsec - err;
}


namespace sfhelp {
	double LepSF()         { return myHelper.LepSF(); }
	double LepSFUp()       { return myHelper.LepSFUp(); }
	double LepSFDown()     { return myHelper.LepSFDown(); }
	double LepSFfastsim()  { return myHelper.LepSFfastsim(); }
	double LepSFfastsimUp(){ return myHelper.LepSFfastsimUp(); }
	double LepSFfastsimDown() { return myHelper.LepSFfastsimDown(); }
	double BtagSF()        { return myHelper.BtagSF(); }
	double BtagHeavyUp()   { return myHelper.BtagHeavyUp(); }
	double BtagHeavyDown() { return myHelper.BtagHeavyDown(); }
	double BtagLightUp()   { return myHelper.BtagLightUp(); }
	double BtagLightDown() { return myHelper.BtagLightDown(); }
	double BtagFSUp()      { return myHelper.BtagFSUp(); }
	double BtagFSDown()    { return myHelper.BtagFSDown(); }
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
	double MetResSF_corr() { return myHelper.MetResSF_corr(); }
	double MetResCorrectionCorridor() { return myHelper.MetResCorrectionCorridor(); }
	double MetResUp()      { return myHelper.MetResUp(); }
	double MetResDown()    { return myHelper.MetResDown(); }
	double TopSystPtSF()   { return myHelper.TopSystPtSF(); }
	double TopSystPtUp()   { return myHelper.TopSystPtUp(); }
	double TopSystPtDown() { return myHelper.TopSystPtDown(); }
	double Contam1lwUp()   { return myHelper.Contam1lwUp(); }
	double Contam1lwDown() { return myHelper.Contam1lwDown(); }
	double ISRnJetsSF()    { return myHelper.ISRnJetsSF(); }
	double ISRnJetsUp()    { return myHelper.ISRnJetsUp(); }
	double ISRnJetsDown()  { return myHelper.ISRnJetsDown(); }
	double LumiUp()        { return myHelper.LumiUp(); }
	double PDFUp()         { return myHelper.PDFUp(); }
	double PDFDown()       { return myHelper.PDFDown(); }
	double WhfXsecUp()     { return myHelper.WhfXsecUp(); }
	double WhfXsecDown()   { return myHelper.WhfXsecDown(); }
	double StopXsecUp()    { return myHelper.StopXsecUp(); }
	double StopXsecDown()  { return myHelper.StopXsecDown(); }
}
