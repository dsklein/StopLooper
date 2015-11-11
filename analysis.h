#ifndef ANALYSIS_H
#define ANALYSIS_H


#include <string>
#include <vector>
#include <iostream>

#include "TString.h"

#include "sample.h"



class analysis{

public:
  void AddSample( sample* newSample );
  void AddSigRegs( std::vector<TString> regions );
  std::vector<short int> GetBkgColors();
  std::vector<TString> GetBkgNamesStorage();
  std::vector<std::string> GetBkgNamesTable();
  std::vector<std::string> GetBkgNamesLegend();
  std::vector<sample*> GetBkgs();
  std::vector<short int> GetSignalColors();
  std::vector<TString> GetSignalNamesStorage();
  std::vector<std::string> GetSignalNamesTable();
  std::vector<std::string> GetSignalNamesLegend();
  std::vector<sample*> GetSignals();
  std::vector<short int> GetColors();
  std::vector<std::vector<TString> > GetSigRegions();
  std::vector<TString> GetSigRegionsAll();
  sample* GetSample( std::string name );
  const int GetNsignals();
  const int GetNbkgs();

private:
  std::vector<sample*> backgrounds;
  std::vector<sample*> signals;
  std::vector<std::vector<TString> > sigRegions;

};

#endif
