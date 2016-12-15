#include "contextVars.h"


// Class "contextvars"
// Provides a set of functions that allow you to set whether variables such as MET, MT2W, etc. should be calculated
// with 2nd lepton added to MET ("_rl") or not, and whether JES should be varied up, down, or not at all.


contextVars myContext; // Declare the extern variable


// Constructor
contextVars::contextVars() {
	useRl = false;
	jesDirection = kNominal;
}

// Other functions
void contextVars::SetJesDir( jesDir direction ) { jesDirection = direction; }
void contextVars::SetUseRl( bool add2ndLep ) { useRl = add2ndLep; }

contextVars::jesDir contextVars::GetJesDir() { return jesDirection; }
bool   contextVars::GetUseRl()  { return useRl; }


const float &contextVars::Met() {
	if(       useRl && jesDirection==kUp      ) return tas::pfmet_rl_jup();
	else if(  useRl && jesDirection==kDown    ) return tas::pfmet_rl_jdown();
	else if(  useRl && jesDirection==kNominal ) return tas::pfmet_rl();
	else if( !useRl && jesDirection==kUp      ) return tas::pfmet_jup();
	else if( !useRl && jesDirection==kDown    ) return tas::pfmet_jdown();
	return tas::pfmet();
}

const float &contextVars::MetPhi() {
	if(       useRl && jesDirection==kUp      ) return tas::pfmet_phi_rl_jup();
	else if(  useRl && jesDirection==kDown    ) return tas::pfmet_phi_rl_jdown();
	else if(  useRl && jesDirection==kNominal ) return tas::pfmet_phi_rl();
	else if( !useRl && jesDirection==kUp      ) return tas::pfmet_phi_jup();
	else if( !useRl && jesDirection==kDown    ) return tas::pfmet_phi_jdown();
	return tas::pfmet_phi();
}

const float &contextVars::MT2W() {
	if(       useRl && jesDirection==kUp      ) return tas::MT2W_rl_jup();
	else if(  useRl && jesDirection==kDown    ) return tas::MT2W_rl_jdown();
	else if(  useRl && jesDirection==kNominal ) return tas::MT2W_rl();
	else if( !useRl && jesDirection==kUp      ) return tas::MT2W_jup();
	else if( !useRl && jesDirection==kDown    ) return tas::MT2W_jdown();
	return tas::MT2W();
}

const float &contextVars::Mindphi_met_j1_j2() {
	if(       useRl && jesDirection==kUp      ) return tas::mindphi_met_j1_j2_rl_jup();
	else if(  useRl && jesDirection==kDown    ) return tas::mindphi_met_j1_j2_rl_jdown();
	else if(  useRl && jesDirection==kNominal ) return tas::mindphi_met_j1_j2_rl();
	else if( !useRl && jesDirection==kUp      ) return tas::mindphi_met_j1_j2_jup();
	else if( !useRl && jesDirection==kDown    ) return tas::mindphi_met_j1_j2_jdown();
	return tas::mindphi_met_j1_j2();
}

const float &contextVars::MT_met_lep() {
	if(       useRl && jesDirection==kUp      ) return tas::mt_met_lep_rl_jup();
	else if(  useRl && jesDirection==kDown    ) return tas::mt_met_lep_rl_jdown();
	else if(  useRl && jesDirection==kNominal ) return tas::mt_met_lep_rl();
	else if( !useRl && jesDirection==kUp      ) return tas::mt_met_lep_jup();
	else if( !useRl && jesDirection==kDown    ) return tas::mt_met_lep_jdown();
	return tas::mt_met_lep();
}

const float &contextVars::TopnessMod() {
	if(       useRl && jesDirection==kUp      ) return tas::topnessMod_rl_jup();
	else if(  useRl && jesDirection==kDown    ) return tas::topnessMod_rl_jdown();
	else if(  useRl && jesDirection==kNominal ) return tas::topnessMod_rl();
	else if( !useRl && jesDirection==kUp      ) return tas::topnessMod_jup();
	else if( !useRl && jesDirection==kDown    ) return tas::topnessMod_jdown();
	return tas::topnessMod();
}

/////////////////////////////////////////////////////////////////////////////////


const bool &contextVars::filt_fastsimjets() {
	if(      jesDirection==kUp   ) return tas::filt_fastsimjets_jup();
	else if( jesDirection==kDown ) return tas::filt_fastsimjets_jdown();
	return tas::filt_fastsimjets();
}

const int &contextVars::ngoodjets() {
	if(      jesDirection==kUp   ) return tas::jup_ngoodjets();
	else if( jesDirection==kDown ) return tas::jdown_ngoodjets();
	return tas::ngoodjets();
}

const int &contextVars::ngoodbtags() {
	if(      jesDirection==kUp   ) return tas::jup_ngoodbtags();
	else if( jesDirection==kDown ) return tas::jdown_ngoodbtags();
	return tas::ngoodbtags();
}

const float &contextVars::ak4_HT() {
	if(      jesDirection==kUp   ) return tas::jup_ak4_HT();
	else if( jesDirection==kDown ) return tas::jdown_ak4_HT();
	return tas::ak4_HT();
}

const float &contextVars::ak4_htratiom() {
	if(      jesDirection==kUp   ) return tas::jup_ak4_htratiom();
	else if( jesDirection==kDown ) return tas::jdown_ak4_htratiom();
	return tas::ak4_htratiom();
}


/////////////////////////////////////////////////////////////////////////////////

const std::vector<float> &contextVars::dphi_ak4pfjet_met() {
	if(      jesDirection==kUp   ) return tas::jup_dphi_ak4pfjet_met();
	else if( jesDirection==kDown ) return tas::jdown_dphi_ak4pfjet_met();
	return tas::dphi_ak4pfjet_met();
}

const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &contextVars::ak4pfjets_p4() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_p4();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_p4();
	return tas::ak4pfjets_p4();
}

const std::vector<bool> &contextVars::ak4pfjets_passMEDbtag() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_passMEDbtag();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_passMEDbtag();
	return tas::ak4pfjets_passMEDbtag();
}

const std::vector<float> &contextVars::ak4pfjets_CSV() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_CSV();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_CSV();
	return tas::ak4pfjets_CSV();
}

const std::vector<float> &contextVars::ak4pfjets_mva() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_mva();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_mva();
	return tas::ak4pfjets_mva();
}

const std::vector<int> &contextVars::ak4pfjets_parton_flavor() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_parton_flavor();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_parton_flavor();
	return tas::ak4pfjets_parton_flavor();
}

const std::vector<int> &contextVars::ak4pfjets_hadron_flavor() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_hadron_flavor();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_hadron_flavor();
	return tas::ak4pfjets_hadron_flavor();
}

const std::vector<bool> &contextVars::ak4pfjets_loose_puid() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_loose_puid();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_loose_puid();
	return tas::ak4pfjets_loose_puid();
}

const std::vector<bool> &contextVars::ak4pfjets_loose_pfid() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_loose_pfid();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_loose_pfid();
	return tas::ak4pfjets_loose_pfid();
}

const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &contextVars::ak4pfjets_leadMEDbjet_p4() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_leadMEDbjet_p4();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_leadMEDbjet_p4();
	return tas::ak4pfjets_leadMEDbjet_p4();
}

const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &contextVars::ak4pfjets_leadbtag_p4() {
	if(      jesDirection==kUp   ) return tas::jup_ak4pfjets_leadbtag_p4();
	else if( jesDirection==kDown ) return tas::jdown_ak4pfjets_leadbtag_p4();
	return tas::ak4pfjets_leadbtag_p4();
}

const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &contextVars::ak4genjets_p4() {
	if(      jesDirection==kUp   ) return tas::jup_ak4genjets_p4();
	else if( jesDirection==kDown ) return tas::jdown_ak4genjets_p4();
	return tas::ak4genjets_p4();
}






namespace context {
	const float &Met() { return myContext.Met(); }
	const float &MetPhi() { return myContext.MetPhi(); }
	const float &MT2W() { return myContext.MT2W(); }
	const float &Mindphi_met_j1_j2() { return myContext.Mindphi_met_j1_j2(); }
	const float &MT_met_lep() { return myContext.MT_met_lep(); }
	const float &TopnessMod() { return myContext.TopnessMod(); }
	const bool &filt_fastsimjets() { return myContext.filt_fastsimjets(); }
	const int &ngoodjets() { return myContext.ngoodjets(); }
	const int &ngoodbtags() { return myContext.ngoodbtags(); }
	const float &ak4_HT() { return myContext.ak4_HT(); }
	const float &ak4_htratiom() { return myContext.ak4_htratiom(); }
	const std::vector<float> &dphi_ak4pfjet_met() { return myContext.dphi_ak4pfjet_met(); }
	const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &ak4pfjets_p4() { return myContext.ak4pfjets_p4(); }
	const std::vector<bool> &ak4pfjets_passMEDbtag() { return myContext.ak4pfjets_passMEDbtag(); }
	const std::vector<float> &ak4pfjets_CSV() { return myContext.ak4pfjets_CSV(); }
	const std::vector<float> &ak4pfjets_mva() { return myContext.ak4pfjets_mva(); }
	const std::vector<int> &ak4pfjets_parton_flavor() { return myContext.ak4pfjets_parton_flavor(); }
	const std::vector<int> &ak4pfjets_hadron_flavor() { return myContext.ak4pfjets_hadron_flavor(); }
	const std::vector<bool> &ak4pfjets_loose_puid() { return myContext.ak4pfjets_loose_puid(); }
	const std::vector<bool> &ak4pfjets_loose_pfid() { return myContext.ak4pfjets_loose_pfid(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ak4pfjets_leadMEDbjet_p4() { return myContext.ak4pfjets_leadMEDbjet_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ak4pfjets_leadbtag_p4() { return myContext.ak4pfjets_leadbtag_p4(); }
	const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &ak4genjets_p4() { return myContext.ak4genjets_p4(); }
}
