#ifndef ANALYSIS_H
#define ANALYSIS_H


#include <string>
#include <vector>
#include <iostream>

#include "TString.h"

#include "sample.h"



class analysis{

public:
  analysis( float lumi, TString fname );

  void AddSample( sample* newSample );
  sample* AddSample( std::string myLabel, std::string niceName );
  sample* AddSample( std::string myLabel, std::string niceName, short int color, sample::sampleType type );
  sample* AddSample( std::string myLabel, std::string tabName, std::string legName );
  sample* AddSample( std::string myLabel, std::string tabName, std::string legName, short int color, sample::sampleType type );
  void AddSigRegs( std::vector<TString> regions );
  std::vector<short int> GetBkgColors();
  std::vector<TString> GetBkgLabels();
  std::vector<std::string> GetBkgNamesTable();
  std::vector<std::string> GetBkgNamesLegend();
  std::vector<sample*> GetBkgs();
  std::vector<short int> GetSignalColors();
  std::vector<TString> GetSignalLabels();
  std::vector<std::string> GetSignalNamesTable();
  std::vector<std::string> GetSignalNamesLegend();
  std::vector<sample*> GetSignals();
  sample* GetData();
  std::vector<short int> GetColors();
  std::vector<std::vector<TString> > GetSigRegions();
  std::vector<TString> GetSigRegionsAll();
  sample* GetSample( std::string name );
  const int GetNsignals();
  const int GetNbkgs();
  bool HasData();
  const float GetLumi();
  const TString GetFileName();

private:
  sample* data;
  std::vector<sample*> backgrounds;
  std::vector<sample*> signals;
  std::vector<std::vector<TString> > sigRegions;
  float luminosity;
  TString plotfilename;
};

#endif
