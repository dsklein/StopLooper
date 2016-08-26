#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H

#include "TString.h"

#include "sfHelper.h"


// Systematic class:
// Holds a reweighting function and a direction (up down, or variation)


class systematic {

public:

  enum direction{ kUp, kDown, kSkipUp, kSkipDown, kVariation };
  systematic( TString systName, direction whichDir, double (*func)() );
  // systematic( TString systName, direction whichDir, double* myWeightVar );
  direction GetDir();
  TString GetName();
  double GetWeight();
  bool IsUp();
  bool IsDown();
  bool IsVariation();
  bool IsSkip();

  static std::vector<TString> GetSysts();
  static std::map<TString,std::vector<direction> > GetVariations();


private:
  TString name;
  direction var_dir;
  double (*reweight_func)();
  // double* weightVar;

  static std::map<TString,std::vector<direction> > variations;

};


#endif
