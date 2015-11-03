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
  std::vector<short int> GetBkgColors();
  std::vector<TString> GetBkgNamesStorage();
  std::vector<std::string> GetBkgNamesTable();
  std::vector<std::string> GetBkgNamesLegend();
  std::vector<short int> GetSignalColors();
  std::vector<TString> GetSignalNamesStorage();
  std::vector<std::string> GetSignalNamesTable();
  std::vector<std::string> GetSignalNamesLegend();
  std::vector<short int> GetColors();
  sample* GetSample( std::string name );
  const int GetNsignals();
  const int GetNbkgs();

private:
  std::vector<sample*> backgrounds;
  std::vector<sample*> signals;

};

#endif
