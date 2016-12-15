#ifndef CONTEXTVARS_H
#define CONTEXTVARS_H


#include "CMS3.h"

#include "Math/LorentzVector.h"


// Class "contextvars"
// Provides a set of functions that allow you to set whether variables such as MET, MT2W, etc. should be calculated
// with 2nd lepton added to MET ("_rl") or not, and whether JES should be varied up, down, or not at all.


class contextVars {


public:
	enum jesDir{ kNominal, kUp, kDown };

	contextVars();

	void SetJesDir( jesDir direction );
	void SetUseRl(  bool   add2ndLep );

	jesDir GetJesDir();
	bool   GetUseRl();

	const float &Met();
	const float &MetPhi();
	const float &MT2W();
	const float &Mindphi_met_j1_j2();
	const float &MT_met_lep();
	const float &TopnessMod();

	const bool &filt_fastsimjets();
	const int &ngoodjets();
	const int &ngoodbtags();
	const float &ak4_HT();
	const float &ak4_htratiom();

	const std::vector<float> &dphi_ak4pfjet_met();
	const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &ak4pfjets_p4();
	const std::vector<bool> &ak4pfjets_passMEDbtag();
	const std::vector<float> &ak4pfjets_CSV();
	const std::vector<float> &ak4pfjets_mva();
	const std::vector<int> &ak4pfjets_parton_flavor();
	const std::vector<int> &ak4pfjets_hadron_flavor();
	const std::vector<bool> &ak4pfjets_loose_puid();
	const std::vector<bool> &ak4pfjets_loose_pfid();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ak4pfjets_leadMEDbjet_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ak4pfjets_leadbtag_p4();
	const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &ak4genjets_p4();


private:
	bool useRl;
	jesDir jesDirection;

};


extern contextVars myContext;


namespace context {
	const float &Met();
	const float &MetPhi();
	const float &MT2W();
	const float &Mindphi_met_j1_j2();
	const float &MT_met_lep();
	const float &TopnessMod();
	const bool &filt_fastsimjets();
	const int &ngoodjets();
	const int &ngoodbtags();
	const float &ak4_HT();
	const float &ak4_htratiom();
	const std::vector<float> &dphi_ak4pfjet_met();
	const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &ak4pfjets_p4();
	const std::vector<bool> &ak4pfjets_passMEDbtag();
	const std::vector<float> &ak4pfjets_CSV();
	const std::vector<float> &ak4pfjets_mva();
	const std::vector<int> &ak4pfjets_parton_flavor();
	const std::vector<int> &ak4pfjets_hadron_flavor();
	const std::vector<bool> &ak4pfjets_loose_puid();
	const std::vector<bool> &ak4pfjets_loose_pfid();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ak4pfjets_leadMEDbjet_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ak4pfjets_leadbtag_p4();
	const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &ak4genjets_p4();
}



#endif
