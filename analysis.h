#ifndef ANALYSIS_H
#define ANALYSIS_H


#include <string>
#include <vector>
#include <iostream>

#include "TString.h"

#include "sample.h"
#include "sigRegion.h"
#include "systematic.h"


class analysis{

public:
  analysis( float lumi, TString fname_plots, TString fname_systs );

  void AddSample( sample* newSample );
  sample* AddSample( std::string myLabel, std::string niceName );
  sample* AddSample( std::string myLabel, std::string niceName, short int color, sample::sampleType type );
  sample* AddSample( std::string myLabel, std::string tabName, std::string legName );
  sample* AddSample( std::string myLabel, std::string tabName, std::string legName, short int color, sample::sampleType type );
  void AddSigRegs( std::vector<sigRegion> regions );
  void AddSystematics( std::vector<systematic*> systs );
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
  std::vector<std::vector<sigRegion> > GetSigRegions();
  std::vector<sigRegion> GetSigRegionsAll();
  std::vector<std::vector<TString> > GetSigRegionLabels();
  std::vector<TString> GetSigRegionLabelsAll();
  std::vector<systematic*> GetSystematics( bool includeSkips );
  std::map<TString,std::vector<TString> > GetSystMap();
  sample* GetSample( std::string name );
  std::vector<sample*> GetAllSamples();
  const int GetNsignals();
  const int GetNbkgs();
  bool HasData();
  const float GetLumi();
  const TString GetPlotFileName();
  const TString GetSystFileName();
  void SetPlotFileName( TString fname );
  void SetSystFileName( TString fname );
  void ResetSigRegions();

private:
  sample* data;
  std::vector<sample*> backgrounds;
  std::vector<sample*> signals;
  std::vector<std::vector<sigRegion> > sigRegions;
  std::vector<systematic*> syst_vars;
  float luminosity;
  TString plotfilename;
  TString systfilename;
};

#endif
