#include "analysis.h"

// Constructor
analysis::analysis( float lumi = 1.0, TString fname_plots = "default.root", TString fname_systs = "dummy.root" )
	: luminosity(lumi),
	  plotfilename(fname_plots),
	  systfilename(fname_systs)
{
	data = NULL;
}

// Copy functions
analysis* analysis::Copy() { return Copy( luminosity, plotfilename, systfilename ); }

analysis* analysis::Copy( float lumi, TString fname_plots, TString fname_systs ) {
	analysis* mycopy = new analysis( lumi, fname_plots, fname_systs );
	mycopy->data = data;
	mycopy->backgrounds = backgrounds;
	mycopy->signals = signals;
	mycopy->sigRegions = sigRegions;
	mycopy->syst_vars = syst_vars;
	return mycopy;
}

// Everything else

void analysis::AddSample( sample* newSample ) {
	if( newSample->IsData() ) data = newSample;
	else if( newSample->IsSignal() ) signals.push_back( newSample );
	else backgrounds.push_back( newSample );
}

sample* analysis::AddSample(std::string myLabel, std::string niceName) {
	sample* myPointer = new sample( myLabel, niceName );
	AddSample( myPointer );
	return myPointer;
}

sample* analysis::AddSample(std::string myLabel, std::string niceName, short int color, sample::sampleType type = sample::kBackground) {
	sample* myPointer = new sample( myLabel, niceName, color, type );
	AddSample( myPointer );
	return myPointer;
}

sample* analysis::AddSample(std::string myLabel, std::string tabName, std::string legName) {
	sample* myPointer = new sample( myLabel, tabName, legName );
	AddSample( myPointer );
	return myPointer;
}

sample* analysis::AddSample(std::string myLabel, std::string tabName, std::string legName, short int color, sample::sampleType type = sample::kBackground) {
	sample* myPointer = new sample( myLabel, tabName, legName, color, type );
	AddSample( myPointer );
	return myPointer;
}

void analysis::AddSigRegs( std::vector<sigRegion> regions ) { sigRegions.push_back(regions); }

void analysis::AddSystematics( std::vector<systematic*> systs ) { syst_vars.insert( syst_vars.end(), systs.begin(), systs.end() ); }

std::vector<short int> analysis::GetBkgColors() {
	std::vector<short int> colors;
	for( sample* mySample : backgrounds ) colors.push_back( mySample->GetColor() );
	return colors;
}

std::vector<TString> analysis::GetBkgLabels() {
	std::vector<TString> labels;
	for( sample* mySample : backgrounds ) labels.push_back( mySample->GetLabel() );
	return labels;
}

std::vector<std::string> analysis::GetBkgNamesTable() {
	std::vector<std::string> tabNames;
	for( sample* mySample : backgrounds ) tabNames.push_back(  mySample->GetTableName().Data() );
	return tabNames;
}

std::vector<std::string> analysis::GetBkgNamesLegend() {
	std::vector<std::string> legNames;
	for( sample* mySample : backgrounds ) legNames.push_back(  mySample->GetLegName().Data() );
	return legNames;
}

std::vector<sample*> analysis::GetBkgs() { return backgrounds; }

std::vector<short int> analysis::GetSignalColors() {
	std::vector<short int> colors;
	for( sample* mySample : signals ) colors.push_back( mySample->GetColor() );
	return colors;
}

std::vector<TString> analysis::GetSignalLabels() {
	std::vector<TString> labels;
	for( sample* mySample : signals ) labels.push_back( mySample->GetLabel() );
	return labels;
}

std::vector<std::string> analysis::GetSignalNamesTable() {
	std::vector<std::string> tabNames;
	for( sample* mySample : signals ) tabNames.push_back(  mySample->GetTableName().Data() );
	return tabNames;
}

std::vector<std::string> analysis::GetSignalNamesLegend() {
	std::vector<std::string> legNames;
	for( sample* mySample : signals ) legNames.push_back(  mySample->GetLegName().Data() );
	return legNames;
}

std::vector<sample*> analysis::GetSignals() { return signals; }

sample* analysis::GetData() {
	if( data ) return data;
	std::cout << "Error: no data sample defined for this analysis!" << std::endl;
	throw(5);
}

std::vector<short int> analysis::GetColors() {
	std::vector<short int> colors;
	for( sample* mySample : backgrounds ) colors.push_back( mySample->GetColor() );
	for( sample* mySample : signals     ) colors.push_back( mySample->GetColor() );
	return colors;
}

sample* analysis::GetSample( std::string name ) {
	for( sample* mySample : backgrounds ) if( mySample->GetLabel() == name ) return mySample;
	for( sample* mySample : signals     ) if( mySample->GetLabel() == name ) return mySample;
	if( data && data->GetLabel() == name ) return data;
	std::cout << "Error! Sample '" << name << "' was not found!" << std::endl;
	throw(5);
}

std::vector<sample*> analysis::GetAllSamples() {
	std::vector<sample*> mySamples;
	if( data ) mySamples.push_back( data );
	mySamples.insert( mySamples.end(), signals.begin(), signals.end() );
	mySamples.insert( mySamples.end(), backgrounds.begin(), backgrounds.end() );
	return mySamples;
}

const int analysis::GetNsignals() { return static_cast<int>(signals.size()); }
const int analysis::GetNbkgs() { return static_cast<int>(backgrounds.size()); }

bool analysis::HasData() {
	if( data ) return true;
	else return false;
}

std::vector<std::vector<sigRegion> > analysis::GetSigRegions() { return sigRegions; }

std::vector<sigRegion> analysis::GetSigRegionsAll() {
	std::vector<sigRegion> output;
	for( std::vector<sigRegion> SRset : sigRegions ) output.insert( output.end(), SRset.begin(), SRset.end() );
	return output;
}

std::vector<std::vector<TString> > analysis::GetSigRegionLabels() {
	std::vector<std::vector<TString> > output;

	for( std::vector<sigRegion> regList : sigRegions ) {
		std::vector<TString> labels;
		for( sigRegion myReg : regList ) labels.push_back( myReg.GetLabel() );
		output.push_back( labels );
	}
	return output;
}

std::vector<TString> analysis::GetSigRegionLabelsAll() {
	std::vector<TString> output;
	for( std::vector<sigRegion> regList : sigRegions ) {
		for( sigRegion myReg : regList ) output.push_back( myReg.GetLabel() );
	}
	return output;
}

std::vector<systematic*> analysis::GetSystematics( bool includeSkips = false ) {
	if( includeSkips ) return syst_vars;
	std::vector<systematic*> syst_list;
	for( systematic* thisvar : syst_vars )  if( !thisvar->IsSkip() ) syst_list.push_back( thisvar );
	return syst_list;
}

std::map<TString,std::vector<TString> > analysis::GetSystMap() {
	std::map<TString,std::vector<TString> > systmap;
	for( systematic* thisVar : syst_vars ) systmap[thisVar->GetName()].push_back(thisVar->GetNameLong());
	return systmap;
}

const float analysis::GetLumi() { return luminosity; }

const TString analysis::GetPlotFileName() { return plotfilename; }

const TString analysis::GetSystFileName() { return systfilename; }

void analysis::SetPlotFileName( TString fname ) { plotfilename = fname; }

void analysis::SetSystFileName( TString fname ) { systfilename = fname; }

void analysis::ResetSigRegions() { sigRegions.clear(); }
