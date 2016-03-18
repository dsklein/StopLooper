#include "analysis.h"

// Constructor
analysis::analysis( float lumi)
  : luminosity(lumi)
{
  data = NULL;
}

// Everything else

void analysis::AddSample( sample* newSample ) {
  if( newSample->IsData() ) data = newSample;
  else if( newSample->IsSignal() ) signals.push_back( newSample );
  else backgrounds.push_back( newSample );
}

void analysis::AddSigRegs( std::vector<TString> regions ) { sigRegions.push_back(regions); }

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

const int analysis::GetNsignals() { return static_cast<int>(signals.size()); }
const int analysis::GetNbkgs() { return static_cast<int>(backgrounds.size()); }

bool analysis::HasData() {
  if( data ) return true;
  else return false;
}

std::vector<std::vector<TString> > analysis::GetSigRegions() { return sigRegions; }

std::vector<TString> analysis::GetSigRegionsAll() {
  std::vector<TString> output;
  for( std::vector<TString> SRset : sigRegions ) output.insert( output.end(), SRset.begin(), SRset.end() );
  return output;
}

const float analysis::GetLumi() { return luminosity; }
