#include "analysis.h"


void analysis::AddSample( sample* newSample ) {
  if( newSample->IsData() ) signals.push_back( newSample );
  else                  backgrounds.push_back( newSample );
}

std::vector<short int> analysis::GetBkgColors() {
  std::vector<short int> colors;
  for( sample* mySample : backgrounds ) colors.push_back( mySample->GetColor() );
  return colors;
}

std::vector<TString> analysis::GetBkgNamesStorage() {
  std::vector<TString> stoNames;
  for( sample* mySample : backgrounds ) stoNames.push_back( mySample->GetIntName() );
  return stoNames;
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

std::vector<short int> analysis::GetSignalColors() {
  std::vector<short int> colors;
  for( sample* mySample : signals ) colors.push_back( mySample->GetColor() );
  return colors;
}

std::vector<TString> analysis::GetSignalNamesStorage() {
  std::vector<TString> stoNames;
  for( sample* mySample : signals ) stoNames.push_back( mySample->GetIntName() );
  return stoNames;
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

std::vector<short int> analysis::GetColors() {
  std::vector<short int> colors;
  for( sample* mySample : backgrounds ) colors.push_back( mySample->GetColor() );
  for( sample* mySample : signals     ) colors.push_back( mySample->GetColor() );
  return colors;
}

sample* analysis::GetSample( std::string name ) {
  for( sample* mySample : backgrounds ) if( mySample->GetIntName() == name ) return mySample;
  for( sample* mySample : signals     ) if( mySample->GetIntName() == name ) return mySample;
  std::cout << "Error! Sample '" << name << "' was not found!" << std::endl;
  throw(5);
}

const int analysis::GetNsignals() { return static_cast<int>(signals.size()); };
const int analysis::GetNbkgs() { return static_cast<int>(backgrounds.size()); };
